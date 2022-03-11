import astropy.units as un
import astropy.coordinates as co
import gammapy.analysis as ga
import gammapy.data as gd
import gammapy.datasets as gds
import gammapy.estimators as ge
import gammapy.makers as gm
import gammapy.maps as gmap
import gammapy.modeling as gmo
import numpy as np
import regions as reg
import os
import shutil
import yaml


def print_diagnostics(run_table):
   tab = run_table.copy()
   tab["excess"].format = "8.5f"
   tab["sqrt_ts"].format = "8.5f"
   tab["alpha"].format = "8.5f"
   tab[["name","counts","counts_off","excess","sqrt_ts","alpha"]].pprint()


def get_exported_observations(src_dir,include_only=None,official_prod=None):
    obs_list = []
    obs_ids = []
    DSs = []
    num_obs = {}
    if not official_prod:
       data_dirs = [loc for loc in os.listdir(src_dir) if "out" in loc]
    else:
       data_dirs = [official_prod]

    for dr in data_dirs:
        ds = gd.DataStore.from_dir(f"{src_dir}/{dr}")
        DSs.append(ds)
        obs_ids.append(ds.obs_table["OBS_ID"].data)

    if include_only:
        print("Only getting: ",include_only)
        for idx,lst in enumerate(obs_ids):
            obs_ids[idx] = [idx for idx in lst if idx in include_only]

    for idx,ds in enumerate(DSs):
        for run in obs_ids[idx]:            
            obs_list.append(ds.obs(run))
        num_obs[data_dirs[idx][4:]] =  len(obs_ids[idx])
    
    for key in num_obs.keys():
       print(f"Found {num_obs[key]} runs in runlist {key}")
    return gd.Observations(obs_list),obs_ids

def map_axis_from_config(axes_config,name):
    return gmap.MapAxis.from_energy_bounds(axes_config.min,
                                           axes_config.max,
                                           nbin=axes_config.nbins,
                                           unit=axes_config.min.unit,
                                           name=name)

def make_exclude_map_from_skypos(source_pos,conf):
    src_excl = reg.CircleSkyRegion(
        center=source_pos,
        radius= (conf["exclude_radius"])*un.deg
    )
    exclude_map = gmap.Map.create(
        npix=(170,170),
        binsz=0.01,
        skydir=source_pos,
        frame="icrs")
    exclude_map = exclude_map.geom.to_image().region_mask([src_excl],inside=False)
    return exclude_map

def make_Analysis_config_from_yaml(config_file):
    """load my custom config file"""
    if isinstance(config_file, dict):
        d = config_file
    else:
        with open(config_file) as f:
            d = yaml.safe_load(f)
    conf = d["mandatory"]
    src_pos = co.SkyCoord(conf["source_pos"]["ra"],
                             conf["source_pos"]["dec"],
                             unit="deg",
                             frame=conf["source_pos"]["frame"])
    if "optional" in d:
      conf["optional"] = d["optional"]
    else:
      conf["optional"] = {}
    gp_conf = d.get("gammapy_arguments")
    if not gp_conf:
        gp_conf = {}
    if conf["optional"].get("energy_axis",None):
        e_ax = conf["optional"]["energy_axis"]
    else:
        e_ax = {"min":0.1,"max":100, "nbins": 72}
    if conf["optional"].get("energy_axis_true",None):
        et_ax = conf["optional"]["energy_axis_true"]
    else:
        et_ax = {"min":0.1,"max":100, "nbins": 72}

    ana_conf = ga.AnalysisConfig(**gp_conf)

    ana_conf.observations.datastore = conf["data_directory"]
    ana_conf.observations.required_irf = ["aeff","edisp","psf"]

    ana_conf.datasets.containment_correction = conf["containment_correction"]
    ana_conf.datasets.map_selection = ["counts","exposure","edisp"]
    ana_conf.datasets.type = "1d"
    ana_conf.datasets.safe_mask.methods = ["edisp-bias"]
    ana_conf.datasets.safe_mask.parameters = {
          "bias_percent":10.}
    ana_conf.datasets.stack = True

    ana_conf.datasets.geom.wcs.skydir = {
        "frame":"galactic",
        "lat":src_pos.galactic.b,
        "lon":src_pos.galactic.l,
    }
    ana_conf.datasets.on_region = {
        "frame":"galactic",
        "lat":src_pos.galactic.b,
        "lon":src_pos.galactic.l,
        "radius": conf["on_radius"]*un.deg
    }
    ana_conf.datasets.geom.axes.energy = {
        "min":e_ax["min"]*un.TeV,
        "max":e_ax["max"]*un.TeV,
        "nbins":e_ax["nbins"]
    }
    ana_conf.datasets.geom.axes.energy_true = {
        "min":et_ax["min"]*un.TeV,
        "max":et_ax["max"]*un.TeV,
        "nbins":et_ax["nbins"]
    }
    ana_conf.flux_points.source = conf["source"]
    ana_conf.flux_points.energy = {
        "min": conf["fit_min"]*un.TeV,
        "max": conf["fit_max"]*un.TeV,
        "nbins":conf["fit_nbins"],
    }
    ana_conf.flux_points.parameters = {}
    ana_conf.fit.fit_range = {
        "min": conf["fit_min"]*un.TeV,
        "max": conf["fit_max"]*un.TeV
    }

    return ana_conf, src_pos, conf

def check_paths(conf):
    if not os.path.isdir(conf["data_directory"]):
        raise ValueError(f"Could not find directory {conf['data_directory']}")
    if not os.path.isdir(conf["out_path"]):
        os.mkdir(conf["out_path"])

def setup_makers(config,src_pos,extra_conf):
   
    e_reco = map_axis_from_config(config.datasets.geom.axes.energy,"energy")
    e_true = map_axis_from_config(config.datasets.geom.axes.energy_true,"energy_true")

    on_region = reg.CircleSkyRegion(center=src_pos,radius=config.datasets.on_region.radius)

    geom = gmap.RegionGeom.create(region=on_region, axes=[e_reco])
    scaffold = gds.SpectrumDataset.create( geom=geom,
                                      energy_axis_true=e_true)

    exclude_map = make_exclude_map_from_skypos(src_pos,extra_conf)
    bkg_maker = gm.ReflectedRegionsBackgroundMaker(exclusion_mask=exclude_map)

    dm_maker = gm.SpectrumDatasetMaker(selection=config.datasets.map_selection,
                                       containment_correction=config.datasets.containment_correction,
                                       use_region_center=True)

    safe_mask_maker = gm.SafeMaskMaker(methods=config.datasets.safe_mask.methods,
                                       **config.datasets.safe_mask.parameters)

    return bkg_maker,safe_mask_maker,dm_maker,scaffold

def make_fluxpoints(analysis,conf):
    energy_edges = np.logspace(np.log10(analysis.config.fit.fit_range.min.value),
                               np.log10(analysis.config.fit.fit_range.max.value),
                               conf["fit_nbins"]) * un.TeV

    flux_points_est = (ge.FluxPointsEstimator(energy_edges=energy_edges,
                                        source=analysis.config.flux_points.source,
                                        selection_optional="all"
                                        )).run(datasets=analysis.datasets)
    flux_points = gds.FluxPointsDataset(models=analysis.models,data=flux_points_est)
    return flux_points,flux_points_est

def save_stats(dataset,conf):
    run_table = gds.Datasets(dataset).info_table(cumulative=False)
    info_table = gds.Datasets(dataset).info_table(cumulative=True)
    info_table["num runs"] = len(dataset)

    # ### Run by run table
    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_stats_byrun.ecsv'
    print(f"Saving run-by-run stats table at {stats_table_path}")
    run_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)

    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_stats_cumul.ecsv'
    print(f"Saving cumulative stats table at {stats_table_path}")
    info_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)
    return info_table,run_table

def save_fit_result(fit_result,ana_conf,conf,E_thr):

    fit_result_path = f'{conf["out_path"]}/{conf["source"]}_fit_result.ecsv'
    print(f"Saving fit result table at {fit_result_path}")
    fit_result.add_row(
          ["","min_energy",ana_conf.config.fit.fit_range.min,"TeV",0.0,None,None,True,""])
    fit_result.add_row(
          ["","max_energy",ana_conf.config.fit.fit_range.max,"TeV",0.0,None,None,True,""])
    fit_result.add_row(
          ["","treshold_energy",E_thr,"TeV",0.0,None,None,True,""])

    fit_result[["name","value","error","unit"]].write(
          fit_result_path, format="ascii.ecsv",overwrite=True)

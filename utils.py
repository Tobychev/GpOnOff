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
import pathlib as pt

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
       if isinstance(include_only,str):
          print("Only getting runs in: ",include_only)
          included_runs = np.loadtxt(include_only,dtype=np.int64)
          included_runs = included_runs[:,0]

       if isinstance(include_only,list):
           print("Only getting: ",include_only)
           included_runs = include_only

       for idx,lst in enumerate(obs_ids):
          obs_ids[idx] = [idx for idx in lst if idx in included_runs]

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

    conf["optional"] = d.get("optional",{})
    conf["datasets"] = d.get("datasets",{})
    gp_conf = d.get("gammapy_arguments",{})

    e_ax =  conf["optional"].get("energy_axis", {"min":0.25,"max":100, "nbins": 72})
    et_ax =  conf["optional"].get("energy_axis_true", {"min":0.1,"max":120, "nbins": 200})

    safe_mask = {}
    if conf["optional"].get("safe_mask_method",None):
       safe_mask["method"] = conf["optional"]["safe_mask_method"]
       safe_mask["parameters"] = conf["optional"]["safe_mask_parameters"]
    else:
       safe_mask["method"] = "edisp-bias"
       safe_mask["parameters"] = {"bias_percent":10.}

    ana_conf = ga.AnalysisConfig(**gp_conf)

    ana_conf.observations.datastore = conf["data_directory"]
    ana_conf.observations.required_irf = ["aeff","edisp","psf"]

    ana_conf.datasets.containment_correction = conf["containment_correction"]
    ana_conf.datasets.map_selection = ["counts","exposure","edisp"]
    ana_conf.datasets.type = "1d"
    ana_conf.datasets.safe_mask.methods = [safe_mask["method"]]
    ana_conf.datasets.safe_mask.parameters = safe_mask["parameters"]
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
        "nbins":conf["flux_nbins"],
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
    outdir = pt.Path(conf["out_path"])
    outdir.mkdir(parents=True,exist_ok=True)

def setup_makers(config,src_pos,extra_conf):

    e_reco = map_axis_from_config(config.datasets.geom.axes.energy,"energy")
    e_true = map_axis_from_config(config.datasets.geom.axes.energy_true,"energy_true")
    rad_bins = gmap.MapAxis.from_bounds(0,5,nbin=30,unit="deg",name="rad")

    on_region = reg.CircleSkyRegion(center=src_pos,radius=config.datasets.on_region.radius)

    map_geom = gmap.WcsGeom.create(skydir=src_pos,
                                   binz=0.02, width=(5,5), fram="icrs", axes=[e_reco])
    geom = gmap.RegionGeom.create(region=on_region, axes=[e_reco])
    scaffold = gds.SpectrumDataset.create( geom=geom,
                                      energy_axis_true=e_true)


    empty_map = gmap.MapDataset.create(geom=map_geom, energy_axis_true=e_true, name="empty")
    exclude_map = make_exclude_map_from_skypos(src_pos,extra_conf)
    bkg_maker = gm.ReflectedRegionsBackgroundMaker(exclusion_mask=exclude_map)

    spec_maker = gm.SpectrumDatasetMaker(selection=config.datasets.map_selection,
                                       containment_correction=config.datasets.containment_correction,
                                       use_region_center=True)

    safe_mask_maker = gm.SafeMaskMaker(methods=config.datasets.safe_mask.methods,
                                       **config.datasets.safe_mask.parameters)

    return bkg_maker,safe_mask_maker,spec_maker,scaffold


def make_dataset(observations,ana_conf,src_pos,conf):
    # ## Setup makers
    bkg_maker,safe_mask_maker,dm_maker,scaffold = setup_makers(ana_conf,src_pos,conf)

    # ## Reduce data

    full_data = gds.Datasets()
    safe_data = gds.Datasets()

    print("Doing run:")
    for ob in observations:
        dataset = dm_maker.run(
            scaffold.copy(name=str(ob.obs_id)), ob)
        print(ob.obs_id,end="... ",flush=True)
        dataset_on_off = bkg_maker.run(dataset,ob)

        full_data.append(dataset_on_off.copy(name=f"{ob.obs_id}_full"))
        try:
            safe_on_off = safe_mask_maker.run(dataset_on_off,ob)
        except:
            print(f"skipping {ob.obs_id}")
            continue
        safe_data.append(safe_on_off)
    print("done")
    return safe_data,full_data

def make_fluxpoints(analysis,conf):
    print(f"Calculating flux points")
    energy_edges = np.logspace(np.log10(analysis.config.fit.fit_range.min.value),
                               np.log10(analysis.config.fit.fit_range.max.value),
                               conf["flux_nbins"]+1) * un.TeV

    flux_points_est = (ge.FluxPointsEstimator(energy_edges=energy_edges,
                                        source=analysis.config.flux_points.source,
                                        selection_optional="all"
                                        )).run(datasets=analysis.datasets)
    flux_points = gds.FluxPointsDataset(models=analysis.models,data=flux_points_est)
    return flux_points,flux_points_est

def resample_dataset(spectrum,conf,rebin):
    fit_axis = gmap.MapAxis.from_energy_bounds(
          conf["fit_min"], conf["fit_max"],
          nbin=rebin,
          unit="TeV",
          name="energy")
    if isinstance(spectrum,gds.Dataset):
       return spectrum.resample_energy_axis(fit_axis,"energy")
    # Note the plural-s
    elif isinstance(spectrum,gds.Datasets):
       rebinned = gds.Datasets()
       for data in datasets:
          spec = data.copy()
          spec = spec.resample_energy_axis(fit_axis,data.name)
          rebinned.append(spec)
       return rebinned
    else:
       raise NotImplementedError(f"Not able to resample spectrum of type {type(spectrum)}")

def save_stats(dataset,conf,zenith,name):
    run_table = gds.Datasets(dataset).info_table(cumulative=False)
    info_table = gds.Datasets(dataset).info_table(cumulative=True)
    try:
       info_table["num runs"] = len(dataset)
    except TypeError:
       info_table["num runs"] = 1
    run_table.add_column(zenith,name="zenith")
    info_table.add_column(zenith,name="zenith")

    # ### Run by run table
    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_{name}_byrun.ecsv'
    print(f"Saving run-by-run stats table at {stats_table_path}")
    run_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)

    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_{name}_cumul.ecsv'
    print(f"Saving cumulative stats table at {stats_table_path}")
    info_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)
    return info_table,run_table

def save_fit_result(fit_result,ana_conf,conf,E_thr):

    fit_result_path = f'{conf["out_path"]}/{conf["source"]}_fit_result.ecsv'
    print(f"Saving fit result table at {fit_result_path}")
    fit_result.add_row(
          ["","min_energy",ana_conf.config.fit.fit_range.min,"TeV",0.0,None,None,True,"",""])
    fit_result.add_row(
          ["","max_energy",ana_conf.config.fit.fit_range.max,"TeV",0.0,None,None,True,"",""])
    fit_result.add_row(
          ["","threshold_energy",E_thr,"TeV",0.0,None,None,True,"",""])

    fit_result[["name","value","error","unit"]].write(
          fit_result_path, format="ascii.ecsv",overwrite=True)

def quick_print_stats(stats):
    print(stats[[
       "excess",
       "sqrt_ts",
       "ontime",
       "counts_rate",
       "background_rate"]][-1])

def get_parameter_dict(parameters):
   pars = {}
   for itm in parameters.to_dict():
      pars[itm["name"]] = itm
   return pars

def calc_analytical_decorr(analysis):
   """
   Lifted from Manuel Meyer's code

   Compute decorrelation energy analytically after initial fit
   and refit. This only works for a power law. See Eq. 3 in https://arxiv.org/pdf/0911.4252.pdf
   """
   covar = analysis.models.covariance.get_subcovariance(['index', 'amplitude', 'reference']).data
   e_decorr = np.exp(covar[0,1] / covar[0,0] / analysis.models.parameters["amplitude"].value)
   e_decorr *= analysis.models.parameters["reference"].value
   return e_decorr

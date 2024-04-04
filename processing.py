import yaml
import regions as reg
import astropy.coordinates as co
import astropy.units as un
import gammapy.analysis as ga
import gammapy.datasets as gds
import gammapy.estimators as ges 
import gammapy.makers as gm
import gammapy.maps as gmap

def map_axis_from_config(axes_config, name):
    return gmap.MapAxis.from_energy_bounds(
        axes_config.min,
        axes_config.max,
        nbin=axes_config.nbins,
        unit=axes_config.min.unit,
        name=name,
    )

def make_exclude_map_from_skypos(source_pos, conf):
    src_excl = reg.CircleSkyRegion(
        center=source_pos, radius=(conf["exclude_radius"]) * un.deg
    )
    exclude_map = gmap.Map.create(
        npix=(170, 170), binsz=0.01, skydir=source_pos, frame="icrs"
    )
    exclude_map = exclude_map.geom.to_image().region_mask([src_excl], inside=False)
    return exclude_map

def maker_common_setup(config, src_pos, extra_conf):
    e_reco = map_axis_from_config(config.datasets.geom.axes.energy, "energy")
    e_true = map_axis_from_config(config.datasets.geom.axes.energy_true, "energy_true")

    on_region = reg.CircleSkyRegion(
        center=src_pos, radius=config.datasets.on_region.radius
    )

    geom = gmap.RegionGeom.create(region=on_region, axes=[e_reco])

    exclude_map = make_exclude_map_from_skypos(src_pos, extra_conf)

    safe_mask_maker = gm.SafeMaskMaker(
        methods=config.datasets.safe_mask.methods,
        **config.datasets.safe_mask.parameters,
    )

    map_geom = gmap.WcsGeom.create(
        skydir=src_pos, binsz=0.02, width=(5, 5), frame="icrs", axes=[e_reco]
    )
    empty_map = gds.MapDataset.create(
        geom=map_geom, energy_axis_true=e_true, name="empty"
    )

    return (
        geom,
        e_true,
        exclude_map,
        empty_map,
        safe_mask_maker,
    )

def setup_onoff_makers(config, src_pos, extra_conf):
    geom, e_true, exclude_map, empty_map, safe_mask_maker = maker_common_setup(config, src_pos, extra_conf)

    empty_spec = gds.SpectrumDataset.create(geom=geom, energy_axis_true=e_true)

    bkg_maker = gm.ReflectedRegionsBackgroundMaker(exclusion_mask=exclude_map)

    spec_maker = gm.SpectrumDatasetMaker(
        selection=config.datasets.map_selection,
        containment_correction=config.datasets.containment_correction,
        use_region_center=extra_conf["optional"].get("use_region_center",False),
    )

    map_maker = gm.MapDatasetMaker(selection=["counts", "exposure", "psf"])

    return (
        bkg_maker,
        safe_mask_maker,
        spec_maker,
        empty_spec,
        map_maker,
        empty_map,
    )

def setup_no_bkg_ring_makers(config, src_pos, extra_conf):
    geom, e_true, exclude_map, empty_map, safe_mask_maker = maker_common_setup(config, src_pos, extra_conf)

    inner_r = config.datasets.on_region.radius+0.4*un.rad


    bkg_maker = gm.RingBackgroundMaker(exclusion_mask=exclude_map, 
            r_in=inner_r, width="0.3 deg")


    map_maker = gm.MapDatasetMaker(selection=["counts", "exposure", "edisp", "psf"])

    excess_maker = ges.ExcessMapEstimator(config.datasets.on_region.radius, selection_optional="")

    breakpoint()
    empty_map = gds.MapDataset.from_geoms(
        geom=empty_map.geoms, name="empty"
    )

    return (
        bkg_maker,
        safe_mask_maker,
        map_maker,
        excess_maker,
        empty_map,
    )

def make_Analysis_config_from_yaml(config_file):
    """load my custom config file"""
    if isinstance(config_file, dict):
        d = config_file
    else:
        with open(config_file) as f:
            d = yaml.safe_load(f)
    conf = d["mandatory"]
    src_pos = co.SkyCoord(
        conf["source_pos"]["ra"],
        conf["source_pos"]["dec"],
        unit="deg",
        frame=conf["source_pos"]["frame"],
    )

    conf["optional"] = d.get("optional", {})
    conf["datasets"] = d.get("datasets", {})
    gp_conf = d.get("gammapy_arguments", {})

    e_ax = conf["optional"].get("energy_axis", {"min": 0.25, "max": 100, "nbins": 72})
    et_ax = conf["optional"].get(
        "energy_axis_true", {"min": 0.1, "max": 120, "nbins": 200}
    )

    safe_mask = {}
    if conf["optional"].get("safe_mask_method", None):
        safe_mask["method"] = conf["optional"]["safe_mask_method"]
        safe_mask["parameters"] = conf["optional"]["safe_mask_parameters"]
    else:
        safe_mask["method"] = "edisp-bias"
        safe_mask["parameters"] = {"bias_percent": 10.0}

    ana_conf = ga.AnalysisConfig(**gp_conf)

    ana_conf.observations.datastore = conf["data_directory"]
    ana_conf.observations.required_irf = ["aeff", "edisp", "psf"]

    ana_conf.datasets.containment_correction = conf["containment_correction"]
    ana_conf.datasets.map_selection = ["counts", "exposure", "edisp"]
    ana_conf.datasets.type = "1d"
    ana_conf.datasets.safe_mask.methods = [safe_mask["method"]]
    ana_conf.datasets.safe_mask.parameters = safe_mask["parameters"]
    ana_conf.datasets.stack = True

    ana_conf.datasets.geom.wcs.skydir = {
        "frame": "galactic",
        "lat": src_pos.galactic.b,
        "lon": src_pos.galactic.l,
    }
    ana_conf.datasets.on_region = {
        "frame": "galactic",
        "lat": src_pos.galactic.b,
        "lon": src_pos.galactic.l,
        "radius": conf["on_radius"] * un.deg,
    }
    ana_conf.datasets.geom.axes.energy = {
        "min": e_ax["min"] * un.TeV,
        "max": e_ax["max"] * un.TeV,
        "nbins": e_ax["nbins"],
    }
    ana_conf.datasets.geom.axes.energy_true = {
        "min": et_ax["min"] * un.TeV,
        "max": et_ax["max"] * un.TeV,
        "nbins": et_ax["nbins"],
    }
    ana_conf.flux_points.source = conf["source"]
    ana_conf.flux_points.energy = {
        "min": conf["fit_min"] * un.TeV,
        "max": conf["fit_max"] * un.TeV,
        "nbins": conf["flux_nbins"],
    }
    ana_conf.flux_points.parameters = {}
    ana_conf.fit.fit_range = {
        "min": conf["fit_min"] * un.TeV,
        "max": conf["fit_max"] * un.TeV,
    }

    return ana_conf, src_pos, conf

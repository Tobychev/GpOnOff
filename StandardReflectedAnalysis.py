#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import astropy.table as tab
import gammapy.analysis as ga
import gammapy.data as gda
import gammapy.datasets as gds
import matplotlib.pyplot as pl
import numpy as np

import utils as ut
import visualisation as vis
HESS1_END = 83490
HESS2_END = 123799

ZEN_BINS = [
    0.0,
    10.0,
    20.0,
    30.0,
    40.0,
    45.0,
    50.0,
    55.0,
    60.0,
    63.0,
    65.0,
]


def split_runlist(runlist):
    with open(runlist, "r") as fil:
        runlines = fil.readlines()
    hess1 = []
    hess2 = []
    hess1u = []
    for runline in runlines:
        run = int(runline.split()[0])
        if run <= HESS1_END:
            hess1.append(run)
        elif run <= HESS2_END:
            hess2.append(run)
        else:
            hess1u.append(run)

    return {"hess1": hess1, "hess2": hess2, "hess1u": hess1u}


def flatten_obs_dict(obs_dict):
    obs_lst = []
    names = []
    props = []
    for key in obs_dict:
        obs, prop = obs_dict[key]
        for itm in obs:
            ob, zen_bin = itm
            obs_lst.append(ob)
            names.append(f"{key}_zbin{zen_bin}")
        props.append(prop)

    prop = tab.vstack(props, join_type="exact")
    return obs_lst, prop, names


def make_dataset(observations, name, ana_conf, src_pos, conf):
    # ## Setup makers
    (
        bkg_maker,
        safe_mask_maker,
        spec_maker,
        empty_spec,
        map_maker,
        empty_map,
    ) = ut.setup_makers(ana_conf, src_pos, conf)

    # ## Reduce data

    full_data = gds.Datasets()
    safe_data = gds.Datasets()
    maps = empty_map.copy()

    print("Doing run:")
    for ob in observations:
        dataset = spec_maker.run(empty_spec.copy(name=f"{name}_{ob.obs_id}"), ob)
        print(ob.obs_id, end="... ", flush=True)
        dataset_on_off = bkg_maker.run(dataset, ob)

        full_data.append(dataset_on_off.copy(name=f"{name}_{ob.obs_id}_full"))
        try:
            safe_on_off = safe_mask_maker.run(dataset_on_off, ob)
        except:
            print(f"skipping {ob.obs_id}")
            continue
        safe_data.append(safe_on_off)

        dataset = map_maker.run(empty_map, ob)
        maps.stack(dataset)

    print("done")
    return safe_data, full_data, maps


def get_listed_observations(data_dir, run_list):
    store = gda.DataStore.from_dir(data_dir)
    obs_lists = []

    tel_ids = tab.Table([run_list], names=["OBS_ID"])
    tel_ids.add_index("OBS_ID")

    props = tab.join(store.obs_table, tel_ids, join_type="inner")
    props["ZEN_BIN"] = np.digitize(props["ZEN_PNT"], ZEN_BINS)

    for group in props.group_by("ZEN_BIN").groups:
        obs_list = []
        for obs_id in group["OBS_ID"]:
            #obs_list.append(store.obs(obs_id, required_irf="point-like"))
            obs_list.append(store.obs(obs_id))
        obs_lists.append((obs_list, group["ZEN_BIN"][0]))

    missing = set(tel_ids["OBS_ID"]) - set(props["OBS_ID"])
    if len(missing) > 0:
        print("WARNING: Some runs were not found:")
        print(sorted(missing))

    return obs_lists, props

def get_all_exported(root_dir):
    data_dirs = list(root_dir.glob("out_*"))

    obs = {}
    for dr in data_dirs:
         ds = gda.DataStore.from_dir(dr)
         runs = ds.obs_table["OBS_ID"].data
         key = dr.name.replace("out_","").replace("_runs","")
         obs[key] = get_listed_observations(dr, runs)
    return obs

def get_listed_hap_observations(runlist, cut_config, prod="fits_prod05"):
    root_dir = pt.Path("/home/hfm/hess/fits/hap-hd/" + prod)

    lists = split_runlist(runlist)

    obs = {}
    for key in lists:
        if len(lists[key]) > 0:
            data_dir = root_dir / key / cut_config
            obs[key] = get_listed_observations(data_dir, lists[key])

    return obs


def get_listed_hd_fr_observations(runlist, cut_config, prod="Prod23_Calib0834"):
    lists = split_runlist(runlist)
    obs = {}
    runs = []
    for key in lists:
        runs += lists[key]

    data_dir = pt.Path("/home/hfm/hess/fits/hap-fr/") / prod / cut_config
    obs["hap-fr"] = get_listed_observations(data_dir, runs)
    return obs


def save_stats(datasets, ids, out_loc, conf, name):
    run_table = datasets.info_table(cumulative=False)
    run_table["name"] = ids
    info_table = datasets.info_table(cumulative=True)
    info_table["name"] = ids
    try:
        info_table["num runs"] = len(datasets)
    except TypeError:
        info_table["num runs"] = 1

    # ### Run by run table
    stats_table_path = out_loc / f'{conf["source"]}_{name}_byrun.ecsv'
    print(f"Saving run-by-run stats table at {stats_table_path}")
    run_table
    run_table.write(stats_table_path, format="ascii.ecsv", overwrite=True)

    stats_table_path = out_loc / f'{conf["source"]}_{name}_cumul.ecsv'
    print(f"Saving cumulative stats table at {stats_table_path}")
    info_table.write(stats_table_path, format="ascii.ecsv", overwrite=True)
    return info_table, run_table


# Start of MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to Run On Off analysis on gammapy formatted data"
    )
    parser.add_argument(
        "config", type=str, help="Config file for the analysis to execute"
    )
    args = parser.parse_args()

    # # Config
    src_ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)
    rebin = conf["optional"].get("fit_energy_bins", None)
    root_dir = pt.Path(conf["data_directory"]) / conf["optional"].get("production","")
    out_dir = pt.Path(conf["out_path"])

    if "hap-hd" in conf["data_directory"]:
        obs_dict = get_listed_hap_observations(
            conf["optional"]["runlist"], conf["optional"]["cut_conf"]
        )

    if "hap-fr" in conf["data_directory"]:
        obs_dict = get_listed_hd_fr_observations(
            conf["optional"]["runlist"], conf["optional"]["cut_conf"]
        )

    if "fits_export" in conf["data_directory"]:
        obs_dict = get_all_exported(root_dir)

    # ## Reduce data
    
    full_data = gds.Datasets()
    safe_data = gds.Datasets()
    map_data = gds.Datasets()
    stat_data = gds.Datasets()

    obs_list, prop, names = flatten_obs_dict(obs_dict)

    obs_ids = []
    for lst, name in zip(obs_list, names):
        safe_set, full_set, maps = make_dataset(lst, name, src_ana_conf, src_pos, conf)
        obid =[itm.obs_id for itm in lst]
        for idx,itm in zip(obid,safe_set):
            stat_data.append(itm)
            obs_ids.append(idx)
        safe_data.append(safe_set.stack_reduce(name=f"{name}_safe"))
        full_data.append(full_set.stack_reduce(name=f"{name}_full"))
        map_data.append(maps.to_image(name=f"{name}_maps"))


    out_loc = out_dir / "reduced"
    out_loc.mkdir(parents=True, exist_ok=True)

    save_stats(stat_data, obs_ids, out_dir, conf, "safe_range")
    #save_stats(full_data, out_loc, conf, "full_range")

    safe_data.write(out_loc / f"{conf['source']}_safe_dataset.yaml", overwrite=True)
    full_data.write(out_loc / f"{conf['source']}_full_dataset.yaml", overwrite=True)
    map_data.write(out_loc / f"{conf['source']}_map_dataset.yaml", overwrite=True)

    safe_tot = safe_data.stack_reduce(name=f"{conf['source']}_safe_total")
    full_tot = full_data.stack_reduce(name=f"{conf['source']}_full_total")

    tab.Table([full_tot.info_dict()]).write(
        out_dir / f"{conf['source']}_full_stat.ecsv", overwrite=True
    )
    tab.Table([safe_tot.info_dict()]).write(
        out_dir / f"{conf['source']}_safe_stat.ecsv", overwrite=True
    )
    prop.write(
        out_dir / f"{conf['source']}_obs_prop.ecsv", format="ascii.ecsv", overwrite=True
    )

    vis.plot_source_stat(stat_data.info_table(cumulative=True),path=conf["out_path"],prefix=f"{conf['source']}_safe")


    for sky in map_data:
        fig = pl.figure()
        mid = sky.geoms["geom"].center_skydir      
        ax = sky.cutout(mid,width=3.).counts.smooth("0.05 deg").plot(add_cbar=True)
        ax.scatter(mid.ra.deg,mid.dec.deg,marker="o",edgecolor="white",transform=ax.get_transform("icrs"))
        pl.savefig(out_dir / f"{conf['source']}_{sky.name}_counts.png")
        pl.close()

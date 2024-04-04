import os

import astropy.table as tab
import gammapy.data as gda
import numpy as np
import pathlib as pt

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

def check_paths(conf):
    if not os.path.isdir(conf["data_directory"]):
        raise ValueError(f"Could not find directory {conf['data_directory']}")
    outdir = pt.Path(conf["out_path"])
    outdir.mkdir(parents=True, exist_ok=True)

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

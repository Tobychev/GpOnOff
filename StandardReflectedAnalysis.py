#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import astropy.table as tab
import gammapy.analysis as ga
import gammapy.data as gda
import gammapy.datasets as gds
import numpy as np

import utils as ut
import visualisation as vis

HESS1_END = 83490
HESS2_END = 123799

ZEN_BINS = [0., 10., 20., 30., 40., 45., 50., 55., 60., 63., 65.,]

def split_runlist(runlist):
   with open(runlist,"r") as fil:
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

   return {"hess1":hess1,"hess2":hess2,"hess1u":hess1u}

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

   prop = tab.vstack(props,join_type="exact")
   return obs, prop,names

def make_dataset(observations,name,ana_conf,src_pos,conf):
    # ## Setup makers
    bkg_maker,safe_mask_maker,spec_maker,scaffold = ut.setup_makers(ana_conf,src_pos,conf)

    # ## Reduce data

    full_data = gds.Datasets(name=name)
    safe_data = gds.Datasets(name=namme)

    print("Doing run:")
    for ob in observations:
        dataset = spec_maker.run(
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

def get_listed_observations(data_dir,run_list):

   store = gda.DataStore.from_dir(data_dir)
   obs_lists = []

   tel_ids = tab.Table([run_list],names=["OBS_ID"])
   tel_ids.add_index("OBS_ID")

   props = tab.join(store.obs_table,tel_ids,join_type="inner")
   props["ZEN_BIN"] = np.digitize(props["ZEN_PNT"],ZEN_BINS)

   for group in props.group_by("ZEN_BIN").groups:
      obs_list = []
      for obs_id in group["OBS_ID"]:
            obs_list.append(store.obs(obs_id))
      obs_lists.append((obs_list,group["ZEN_BIN"][0])

   missing = set(tel_ids["OBS_ID"])-set(props["OBS_ID"])
   if len(missing) > 0:
       print("WARNING: Some runs were not found:")
       print(sorted(missing))

   return obs_lists,props


def get_listed_hap_observations(runlist,cut_config,prod="fits_prod05"):
   root_dir = pt.Path("/home/hfm/hess/fits/hap-hd/"+ prod)

   lists = split_runlist(runlist)

   obs = {}
   for key in lists:
      if len(lists[key]) > 0:
         data_dir = root_dir / key / cut_config
         obs[key] = get_listed_observations(data_dir, lists[key])

   return obs


def get_listed_hd_fr_observations(runlist,cut_config,prod="Prod23_Calib0834"):
   lists = split_runlist(runlist)
   obs = {}
   runs = []
   for key in lists:
      runs += lists[key]

   data_dir = pt.Path("/home/hfm/hess/fits/hap-fr/") / prod / cut_config
   obs["hap-fr"] = get_listed_observations(data_dir, runs)
   return obs

def save_stats(datasets,conf,name):
    run_table = datasets.info_table(cumulative=False)
    info_table = datasets.info_table(cumulative=True)
    try:
       info_table["num runs"] = len(datasets)
    except TypeError:
       info_table["num runs"] = 1

    # ### Run by run table
    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_{name}_byrun.ecsv'
    print(f"Saving run-by-run stats table at {stats_table_path}")
    run_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)

    stats_table_path = f'{conf["out_path"]}/{conf["source"]}_{name}_cumul.ecsv'
    print(f"Saving cumulative stats table at {stats_table_path}")
    info_table.write( stats_table_path ,format="ascii.ecsv", overwrite=True)
    return info_table,run_table


# Start of MAIN
if __name__ == "__main__":

   parser = argparse.ArgumentParser(description="Script to Run On Off analysis on gammapy formatted data")
   parser.add_argument("config",
                      type=str,
                      help="Config file for the analysis to execute")
   args = parser.parse_args()

   # # Config
   src_ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
   ut.check_paths(conf)
   rebin = conf["optional"].get("fit_energy_bins",None)
   root_dir = pt.Path(conf["data_directory"]) / conf["optional"]["production"]
   out_dir = pt.Path(conf["out_path"])

   if "hap-hd" in conf["data_directory"]:
      obs_dict = get_listed_hap_observations(conf["optional"]["runlist"],conf["optional"]["cut_conf"])
      datastore = "{root_dir}/{key}/{cut_conf}"

   if "hap-fr" in conf["data_directory"]:
      obs_dict = get_listed_hd_fr_observations(conf["optional"]["runlist"],conf["optional"]["cut_conf"])
      datastore = "{root_dir}/{cut_conf}"

   # ## Reduce data
   full_data = gds.Datasets()
   safe_data = gds.Datasets()

   obs_list, prop,names = flatten_obs_dict(obs_dict)
   for lst,name in zip(obs_list,names):
      safe_set, full_set = make_dataset(lst,name,src_ana_conf,src_pos,conf)
      safe_data.append(safe_set.stack_reduce())
      full_data.append(full_set.stack_reduce())

   save_stats(safe_data,conf,"safe_range")
   save_stats(full_data,conf,"full_range")
   breakpoint()

   out_loc = out_dir / "reduced"
   out_loc.mkdir(parents=True,exist_ok=True)
   safe_data.write(out_loc / f"{conf['source']}_safe_dataset.yaml", overwrite=True)
   safe_data.write(out_loc / f"{conf['source']}_full_dataset.yaml", overwrite=True)

   

   for key in obs_dict:
       obs_objs, obs_tab = obs_dict[key]
       ana_conf = src_ana_conf.copy()

       out_loc = out_dir / key
       out_loc.mkdir(parents=True,exist_ok=True)
       ana_conf.general.datasets_file = out_loc / f"{conf['source']}_{key}_dataset.yaml"

       ana_conf.observations.obs_ids = [itm.obs_id for itm in obs_objs]
       ana_conf.observations.datastore = pt.Path(
          datastore.format(
             root_dir= root_dir,
             key= key,
             cut_conf = conf["optional"]["cut_conf"]))

       ana_conf.observations.required_irf = ["aeff","edisp","psf", "bkg"]

       ana_conf.datasets.map_selection = ["counts","exposure","edisp", "background"]
       ana_conf.datasets.background.method = "reflected"

       print(ana_conf)
       ana = ga.Analysis(ana_conf)
       ana.get_observations()
       ana.get_datasets()

       info = ana.datasets.info_table()


       ana.config.write(out_loc / f"{key}_conf.log", overwrite = True)
       info.write(out_loc / f"{key}_stats.csv", format = "ascii.ecsv", overwrite=True)
       ana.write_datasets()
       print(
         f"Tobs={info['livetime'].to('h')[0]:.1f} Excess={info['excess'].value[0]:.1f} \
               Significance={info['sqrt_ts'][0]:.2f}")



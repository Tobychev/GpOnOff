#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import gammapy.analysis as ga
import gammapy.data as gda
import gammapy.datasets as gds

import utils as ut
import visualisation as vis

HESS1_END = 83490
HESS2_END = 123799

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


def get_listed_observations(data_dir,run_list):

   store = gda.DataStore.from_dir(data_dir)
   obs_list = []
   
   for obs_id in run_list:
      try:
         obs_list.append(store.obs(obs_id))
      except ValueError:
         print(f"could not find {obs_id} in {data_dir}, dropping run")
         run_list.remove(obs_id)
   store.obs_table.add_index("OBS_ID")
   try:
      props = store.obs_table.loc[run_list]
   except KeyError as err:
      print(f"Problem with key: {err}")
      props = {}
   return obs_list,props


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



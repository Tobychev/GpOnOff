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
HESS2_END = 127699

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


def get_listed_hap_observations(runlist,cut_config,prod="fits_prod05"):
   root_dir = pt.Path("/home/hfm/hess/fits/hap-hd/"+ prod)

   lists = split_runlist(runlist)

   obs = {}
   for key in lists:
      if len(lists[key]) > 0:
         data_dir = root_dir / key / cut_config
         store = gda.DataStore.from_dir(data_dir)
         obs_list = []
         for obs_id in lists[key]:
            obs_list.append(store.obs(obs_id))
         store.obs_table.add_index("OBS_ID")
         props = store.obs_table.loc[lists[key]]
         obs[key] = obs_list,props


   return obs

# Start of MAIN
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to Run On Off analysis on gammapy formatted data")
    parser.add_argument("config",
                       type=str,
                       help="Config file for the analysis to execute")
    parser.add_argument("--skip-flux",
          action='store_true',
          default=False,
          help = "Skip calculating fluxes (saves time)")
    parser.add_argument("--debug",
          action='store_true',
          default=False,
          help = "Create extra debug plots")
    args = parser.parse_args()

    # # Config
    src_ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)
    joint_fit = conf["joint_fit"]
    rebin = conf["optional"].get("fit_energy_bins",None)
    root_dir = pt.Path("/home/hfm/hess/fits/hap-hd/") / conf["optional"]["production"] 

    obs_dict = get_listed_hap_observations(conf["optional"]["runlist"],conf["optional"]["cut_conf"])

    for key in obs_dict:
        obs_objs, obs_tab = obs_dict[key]
        ana_conf = src_ana_conf.copy()

        ana_conf.observations.obs_ids = [itm.obs_id for itm in obs_objs]
        ana_conf.observations.datastore = root_dir / key / conf["optional"]["cut_conf"]
        ana_conf.observations.required_irf = ["aeff","edisp","psf", "bkg"]

        ana_conf.datasets.type = "3d"
        ana_conf.datasets.stacked = False

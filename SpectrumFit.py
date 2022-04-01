#!/usr/bin/env python
import argparse
import astropy.units as un
import gammapy.analysis as ga
import gammapy.datasets as gds
import gammapy.maps as gmap
import gammapy.modeling as gmo
import matplotlib.pyplot as pl
import numpy as np
import os
import shutil
import yaml

import utils as ut
import visualisation as vis

def make_model(conf,kind="PL"):
   if kind == "PL":
     model = gmo.models.PowerLawSpectralModel(
          index=2,
          amplitude="1e-12 TeV-1 cm-2 s-1",
          reference=conf["fit_E0"] * un.TeV,
      )
   else:
      raise ValueError(f"Unsuported model type: '{kind}'")
   tot_model = gmo.models.Models([
      gmo.models.SkyModel(
         spectral_model=model,
         name=conf["source"])
      ])
   return tot_model

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Run On Off analysis on gammapy formatted data")
    parser.add_argument("config",
                       type=str,
                       help="Config file for the analysis to execute")
    args = parser.parse_args()

    # # Config
    ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)

    if conf["optional"].get("dataset_file",None):
       print(f'reading data from:\n  {conf["optional"]["dataset_file"]}')
       datasets = gds.Datasets.read(conf["optional"]["dataset_file"])
    else:
       print(f"reading data from:\n  {conf['out_path']}/Dataset/{conf['source']}_dataset.yaml")
       datasets = gds.Datasets.read(
             f"{conf['out_path']}/Dataset/{conf['source']}_dataset.yaml")

   
    # # Fit spectrum 

    ana = ga.Analysis(ana_conf)
    ana.datasets = gds.Datasets(datasets)
    ana.set_models(make_model(conf))

    ana.run_fit()
    ana.fit.covariance(ana.datasets)
    fig, cont = vis.plot_contours(conf,ana,"amplitude","index")
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_Countours.png')

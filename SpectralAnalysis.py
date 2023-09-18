#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import astropy.units as u
import gammapy.analysis as ga
import gammapy.datasets as gds
import gammapy.modeling.models as gmo
import visualisation as vis

import utils as ut

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to Run On Off analysis on gammapy formatted data")
    parser.add_argument("config",
                       type=str,
                       help="Config file for the analysis to execute")
    parser.add_argument("model",
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
    rebin = conf["optional"].get("fit_energy_bins",None)
    out_dir = pt.Path(conf["out_path"])

#    obs_dict = get_listed_hap_observations(conf["optional"]["runlist"],conf["optional"]["cut_conf"])
    
    datasets = gds.Datasets()
    for data in conf["datasets"]:
       dset = gds.Datasets.read(filename=data)
       datasets.append(dset.stack_reduce())

    ana_conf = src_ana_conf.copy()

    ana = ga.Analysis(ana_conf)
    ana.datasets = datasets

    mod = gmo.Models.read(args.model)
    mod.names[0] = conf["source"]
    ana.set_models(mod)
    ana.run_fit()

    Edec = ut.calc_analytical_decorr(ana)
    print(f"Analytical decorrelation energy is {Edec}")
    print("Rerun with this value if desired")

    fit_result = ana.fit_result.parameters.to_table()
    ut.save_fit_result(fit_result,ana,conf,-1.0*u.TeV)

    best_fit = out_dir / f"{conf['source']}_best_fit.yaml"
    ana.models.write(best_fit,overwrite=True)    

    vis.save_fit_residuals(ana.datasets[0],conf,f"PowerLawFit_Residuals")

    if not args.skip_flux:
#       ana.get_flux_points()
       flux_points,flux_points_est = ut.make_fluxpoints(ana,conf)
       vis.save_flux_points(flux_points_est,conf,f"PowerLawFit_Likelihood")
       try:
          vis.save_fluxpoint_fit(flux_points,conf,f"PowerLawFit")
       except ValueError:
          print("Only upper limits found, residuals can't be calculated")

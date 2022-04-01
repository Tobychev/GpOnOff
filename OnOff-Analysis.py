#!/usr/bin/env python
import argparse
import astropy.units as un
import gammapy.analysis as ga
import gammapy.datasets as gds
import gammapy.maps as gmap
import gammapy.modeling as gmo
import numpy as np
import os
import shutil
import yaml

import utils as ut
import visualisation as vis

def make_dataset(observations,ana_conf,src_pos,conf):
    # ## Setup makers
    bkg_maker,safe_mask_maker,dm_maker,scaffold = ut.setup_makers(ana_conf,src_pos,conf)

    # ## Reduce data

    full_data = gds.Datasets()
    safe_data = gds.Datasets()

    print("Doing run:")
    for ob in observations:
        dataset = dm_maker.run(
            scaffold.copy(name=str(ob.obs_id)), ob)
        print(ob.obs_id,end="... ",flush=True)
        dataset_on_off = bkg_maker.run(dataset,ob)
        full_data.append(dataset_on_off.copy())
        try:
            safe_on_off = safe_mask_maker.run(dataset_on_off,ob)
        except:
            print(f"skipping {ob.obs_id}")
            continue
        safe_data.append(safe_on_off)
    print("done")
    return safe_data,full_data


def calc_analytical_decorr(analysis):
   """
   Lifted from Manuel Meyer's code

   Compute decorrelation energy analytically after initial fit
   and refit. This only works for a power law. See Eq. 3 in https://arxiv.org/pdf/0911.4252.pdf
   """
   covar = ana.models.covariance.get_subcovariance(['index', 'amplitude', 'reference']).data
   e_decorr = np.exp(covar[0,1] / covar[0,0] / analysis.models.parameters["amplitude"].value)
   e_decorr *= analysis.models.parameters["reference"].value
   return e_decorr

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
    ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)
    joint_fit = conf["joint_fit"]
    rebin = conf["optional"].get("fit_energy_bins",None)

    # ## Get obs

    observations, obs_ids = ut.get_exported_observations(
             conf["data_directory"],
             conf["optional"].get("include_only",None),
             conf["optional"].get("official_prod",None))
    ana_conf.observations.obs_ids = list(np.concatenate(obs_ids))

    if len(ana_conf.observations.obs_ids) == 0:
      raise ValueError(f"No accepted runs found in {conf['data_directory']}")

    # ## Reduce data

    if conf["optional"].get("dataset_file",None):
       print(f'reading {conf["optional"]["dataset_file"]}')
       datasets = gds.Datasets.read(conf["optional"]["dataset_file"])
       full_data = gds.Datasets.read(
             conf["optional"]["dataset_file"].replace(".yaml","_full.yaml"))
    else:
       datasets,full_data =  make_dataset(observations,ana_conf,src_pos,conf)
       datasets.write(f"{conf['out_path']}/Dataset/{conf['source']}_dataset.yaml",
             overwrite=True)
       full_data.write(f"{conf['out_path']}/Dataset/{conf['source']}_dataset_full.yaml",
             overwrite=True)
    nrun = len(datasets)

    # ## Statistics
    info_table,run_table = ut.save_stats(full_data,conf,"stats")
    ut.quick_print_stats(info_table)
    vis.plot_source_stat(info_table,path=conf["out_path"],prefix=conf["source"])

    # Stack
    stacked_data  = datasets.stack_reduce()
    Ethr_stacked = stacked_data.energy_range_safe[0].data[0][0]

    if args.debug:
       vis.save_exp_corrected_counts(stacked_data,conf,"Initial")
       vis.save_spectrum_diagnostic(stacked_data,conf,"Spectrum_Raw_Binning")


    if joint_fit:
       print("Using joint dataset")
       Ethr = np.zeros(nrun)
       joint_data = datasets
       for idx,data in enumerate(joint_data):
          Ethr[idx] = data.energy_range_safe[0].data[0][0]
       if rebin:
          joint_data = ut.resample_dataset(joint_data,conf,rebin)
          if args.debug:
             rebinned = joint_data.stack_reduce()
             vis.save_spectrum_diagnostic(rebinned,conf,"Spectrum_Rebinned")
             vis.save_exp_corrected_counts(rebinned,conf,"Resampled")

       spectrum = joint_data
    else:
       print("Using stacked dataset")
       if rebin:
          stacked_data = ut.resample_dataset(stacked_data.copy(),conf,rebin)

       spectrum = stacked_data

    # # Fit spectrum 

    ana = ga.Analysis(ana_conf)
    ana.observations = observations
    ana.datasets = gds.Datasets(spectrum)
    pl_model = gmo.models.PowerLawSpectralModel(
        index=2,
        amplitude="1e-12 TeV-1 cm-2 s-1",
        reference=conf["fit_E0"] * un.TeV,
    )
    tot_model = gmo.models.Models([gmo.models.SkyModel(spectral_model=pl_model,name=ana.config.flux_points.source)])

    ana.set_models(tot_model)
    ana.run_fit()
    ana.fit.covariance(ana.datasets)
    Edec = calc_analytical_decorr(ana)
    print(f"Analytical decorrelation energy is {Edec}")
    print("Rerun with this value if desired")

    fit_result = ana.fit_result.parameters.to_table()
    ut.save_fit_result(fit_result,ana,conf,Ethr_stacked)

    if joint_fit:
       fit_type = "Joint"
    else:
       fit_type = "Stacked"
    vis.save_fit_residuals(ana.datasets[0],conf,f"PowerLaw{fit_type}Fit_Residuals")

    if not args.skip_flux:
       flux_points,flux_points_est = ut.make_fluxpoints(ana,conf)
       vis.save_flux_points(flux_points_est,conf,f"PowerLaw{fit_type}Fit_Likelihood")
       vis.save_fluxpoint_fit(flux_points,conf,f"PowerLaw{fit_type}Fit")


    shutil.copy(os.path.abspath(args.config),
                f'{conf["out_path"]}/{conf["source"]}_config.yaml')

    if joint_fit:
       print("Joint fit spectrum stats:")
    else:
       print("Stacked fit spectrum stats:")

    ut.quick_print_stats(
          gds.Datasets(datasets).info_table(
             cumulative=True))

    print(f"Threshold {Ethr_stacked:.4f} \n")
    print(fit_result[[
       "name",
       "value",
       "error",
       "unit",]])
    print()

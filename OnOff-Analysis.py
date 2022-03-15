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

    datasets = gds.Datasets()

    print("Doing run:")
    for ob in observations:
        dataset = dm_maker.run(
            scaffold.copy(name=str(ob.obs_id)), ob)
        print(ob.obs_id,end="... ",flush=True)
        dataset_on_off = bkg_maker.run(dataset,ob)
        try:
            dataset_on_off = safe_mask_maker.run(dataset_on_off,ob)
        except:
            print(f"skipping {ob.obs_id}")
            continue
        datasets.append(dataset_on_off)
    print("done")
    return datasets

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
    else:
       datasets =  make_dataset(observations,ana_conf,src_pos,conf)
       datasets.write(f"{conf['out_path']}/Dataset/{conf['source']}_dataset.yaml",
             overwrite=True)
    nrun = len(datasets)


    # ## Statistics
    info_table,run_table = ut.save_stats(datasets,conf,"stats")
    quick_print_stats(info_table)
    vis.plot_source_stat(info_table,path=conf["out_path"],prefix=conf["source"])

    # Check the binning, only sensible for the stacked dataset
    spectrum  = datasets.stack_reduce()
    vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Raw_Binning")
    Ethr_stacked = spectrum.energy_range_safe[0].data[0][0]
    # Resample to new axis
    rebinned = gds.Datasets()
    if rebin:
       print("Resampling!")
       axis_config = rebin
       fit_axis = gmap.MapAxis.from_energy_bounds(
             rebin["min"], rebin["max"],
             nbin=rebin["nbins"],
             unit="TeV",
             name="energy")
       spectrum = spectrum.resample_energy_axis(fit_axis,"energy")
       vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Rebinned")
       info_table,run_table = ut.save_stats(datasets,conf,"rebinned_stats")

       if joint_fit:
          for data in datasets:
             spec = data.copy()
             spec = spec.resample_energy_axis(fit_axis,data.name)
             rebinned.append(spec)
          print("Using resampled joint dataset")
          spectrum = rebinned

    elif os.path.exists(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png'):
       os.remove(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png')


    if joint_fit and not rebinned:
       print("Using joint dataset")
       spectrum = datasets

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
    print("Fitting to get decorrelaion energy")
    ana.run_fit()
    ana.fit.covariance(ana.datasets)
    Edec = calc_analytical_decorr(ana)


    print("Re-fitting at decorrelaion energy")

    ana.models.parameters["index"].value = 2
    ana.models.parameters["reference"].value = Edec
    ana.models.parameters["amplitude"].value = 1e-12
    ana.models.parameters["reference"].frozen = True
    ana.run_fit()
    ana.fit.covariance(ana.datasets)

    fit_result = ana.fit_result.parameters.to_table()
    ut.save_fit_result(fit_result,ana,conf,Ethr_stacked)

    if joint_fit:
       fit_type = "Joint"
    else:
       fit_type = "Stacked"
    vis.save_fit_residuals(ana.datasets[0],conf,f"PowerLaw{fit_type}Fit_Residuals")

    flux_points,flux_points_est = ut.make_fluxpoints(ana,conf)
    vis.save_flux_points(flux_points_est,conf,f"PowerLaw{fit_type}Fit_Likelihood")
    vis.save_fluxpoint_fit(flux_points,conf,f"PowerLaw{fit_type}Fit")


    shutil.copy(os.path.abspath(args.config),
                f'{conf["out_path"]}/{conf["source"]}_config.yaml')

    if joint_fit:
       print("Joint fit result:")
    else:
       print("Stacked result:")
    quick_print_stats(info_table)

    print(f"Threshold {Ethr_stacked:.4f} \n")
    print(fit_result[[
       "name",
       "value",
       "error",
       "unit",]])
    print()

#!/usr/bin/env python
import argparse
import astropy.units as un
import astropy.coordinates as co
import gammapy.analysis as ga
import gammapy.data as gd
import gammapy.datasets as gds
import gammapy.estimators as ge
import gammapy.makers as gm
import gammapy.maps as gmap
import gammapy.modeling as gmo
import numpy as np
import regions as reg
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to Run On Off analysis on gammapy formatted data")
    parser.add_argument("config",
                       type=str,
                       help="Config file for the analysis to execute")
    parser.add_argument("--save-dataset",
                  action="store_true",
                  help="""Save the dataset made by the On-Off Maker""")
    args = parser.parse_args()

    # # Config
    ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)

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
       datasets = gds.Datasets.read(conf["optional"]["dataset_file"])
    else:
       datasets =  make_dataset(observations,ana_conf,src_pos,conf)
    nrun = len(datasets)

    if args.save_dataset:
       datasets.write(f"{conf['out_path']}/Dataset_{conf['source']}")

    # ## Statistics
    info_table,run_table = ut.save_stats(datasets,conf)
    vis.plot_source_stat(info_table,path=conf["out_path"],prefix=conf["source"])

    # # Fit spectrum stacked

    ana = ga.Analysis(ana_conf)
    ana.observations = observations
    spectrum  = datasets.stack_reduce()
    vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Raw_Binning")

    if conf["optional"].get("fit_energy_axis",None):
       print("Resampling!")
       axis_config = conf["optional"]["fit_energy_axis"]
       fit_axis = gmap.MapAxis.from_energy_bounds(
             axis_config["min"], axis_config["max"],
             nbin=axis_config["nbins"],
             unit="TeV",
             name="energy")
       spectrum = spectrum.resample_energy_axis(fit_axis,"energy")
       vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Rebinned")
       
    ana.datasets = gds.Datasets(spectrum)
    Ethr_stacked = ana.datasets[0].energy_range_safe[0].data[0][0]
    
    pl_model = gmo.models.PowerLawSpectralModel(
        index=2,
        amplitude="1e-12 TeV-1 cm-2 s-1",
        reference=conf["fit_E0"] * un.TeV,
    )
    tot_model = gmo.models.Models([gmo.models.SkyModel(spectral_model=pl_model,name=ana.config.flux_points.source)])
    ana.set_models(tot_model)


    ana.run_fit()
    ana.fit.covariance(ana.datasets)
    fit_result = ana.fit_result.parameters.to_table()
    ut.save_fit_result(fit_result,ana,conf,Ethr_stacked)
    vis.save_fit_residuals(ana.datasets[0],conf,"PowerLawFit_Residuals")

    if conf["optional"].get("fit_energy_axis",None):
       vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Rebinned")
    elif os.path.exists(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png'):
       os.remove(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png')


    flux_points,flux_points_est = ut.make_fluxpoints(ana,conf)
    vis.save_flux_points(flux_points_est,conf,"PowerLawFit_Likelihood")
    vis.save_fluxpoint_fit(flux_points,conf,"PowerLawFit")


    shutil.copy(os.path.abspath(args.config),
                f'{conf["out_path"]}/{conf["source"]}_config.yaml')

    print(info_table[[
       "excess",
       "sqrt_ts",
       "ontime",
       "counts_rate",
       "background_rate"]][-1])
    print()
    print(f"Threshold {Ethr_stacked:.4f} \n")
    print(fit_result[[
       "name",
       "value",
       "error",
       "unit",]])
    print()

"""

    # # Fit spectrum joint

    anaj = ga.Analysis(ana_conf)
    anaj.observations = observations
    spectrum  = datasets.copy
    vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Raw_Binning")

    if conf["optional"].get("fit_energy_axis",None):
       print("Resampling!")
       axis_config = conf["optional"]["fit_energy_axis"]
       fit_axis = gmap.MapAxis.from_energy_bounds(
             axis_config["min"], axis_config["max"],
             nbin=axis_config["nbins"],
             unit="TeV",
             name="energy")
       spectrum = spectrum.resample_energy_axis(fit_axis,"energy")
       vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Rebinned")
       
    ana.datasets = gds.Datasets(spectrum)
    Ethr_stacked = ana.datasets[0].energy_range_safe[0].data[0][0]
    
    pl_model = gmo.models.PowerLawSpectralModel(
        index=2,
        amplitude="1e-12 TeV-1 cm-2 s-1",
        reference=conf["fit_E0"] * un.TeV,
    )
    tot_model = gmo.models.Models([gmo.models.SkyModel(spectral_model=pl_model,name=ana.config.flux_points.source)])
    ana.set_models(tot_model)


    ana.run_fit()
    ana.fit.covariance(ana.datasets)
    fit_result = ana.fit_result.parameters.to_table()
    ut.save_fit_result(fit_result,ana,conf,Ethr_stacked)
    vis.save_fit_residuals(ana.datasets[0],conf,"PowerLawFit_Residuals")

    if conf["optional"].get("fit_energy_axis",None):
       vis.save_spectrum_diagnostic(spectrum,conf,"Spectrum_Rebinned")
    elif os.path.exists(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png'):
       os.remove(f'{conf["out_path"]}/{conf["source"]}_Spectrum_Rebinned.png')

    
    flux_points,flux_points_est = ut.make_fluxpoints(ana,conf)
    vis.save_flux_points(flux_points_est,conf,"PowerLawFit_Likelihood")
    vis.save_fluxpoint_fit(flux_points,conf,"PowerLawFit")
"""

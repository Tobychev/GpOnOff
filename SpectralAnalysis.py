#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import astropy.units as u
import astropy.table as tab
import gammapy.analysis as ga
import gammapy.datasets as gds
import gammapy.modeling.models as gmo
import visualisation as vis

import utils as ut

def make_model(conf,kind="PL"):
   index = conf["optional"].get("start_index",2)
   ampl = conf["optional"].get("start_amplitude",1)
   if kind == "PL":
     model = gmo.PowerLawSpectralModel(
          index=index,
          amplitude=f"{ampl}e-12 TeV-1 cm-2 s-1",
          reference=conf["fit_E0"] * un.TeV,
      )
   else:
      raise ValueError(f"Unsuported model type: '{kind}'")
   tot_model = gmo.Models([
      gmo.SkyModel(
         spectral_model=model,
         name=conf["source"])
      ])
   return tot_model

def set_E_ref_to_pivot(analysis):
    model = analysis.models[0].spectral_model
    if isinstance(model,gmo.PowerLawSpectralModel):
        model.reference.value = ana.models[0].spectral_model.pivot_energy.to_value("TeV")
    else:
        raise ValueError(f"Unsuported model type: '{model}'")
    return analysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to Run On Off analysis on gammapy formatted data"
    )
    parser.add_argument(
        "config", type=str, help="Config file for the analysis to execute"
    )
    parser.add_argument(
        "model", type=str, help="Model file to use when fitting"
    )
    parser.add_argument(
        "--skip-flux",
        action="store_true",
        default=False,
        help="Skip calculating fluxes (saves time)",
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="Create extra debug plots"
    )
    args = parser.parse_args()

    # # Config
    src_ana_conf, src_pos, conf = ut.make_Analysis_config_from_yaml(args.config)
    ut.check_paths(conf)
    rebin = conf["optional"].get("fit_energy_bins", None)
    out_dir = pt.Path(conf["out_path"])

    if conf["optional"].get("contour_points",None):
       npoints = conf["optional"]["contour_points"]
    else:
       npoints = 10

    if len(conf["datasets"]) == 1:
       datasets = gds.Datasets.read(filename=conf["datasets"][0])
    else:
      raise ValueError("Only one dataset file supported now")

    ana_conf = src_ana_conf.copy()

    ana = ga.Analysis(ana_conf)
    ana.datasets = datasets

    if args.model:
      mod = gmo.Models.read(args.model)
      mod.names[0] = conf["source"]
    else:
      mod = make_model(conf)

    ana.set_models(mod)
    ana.run_fit()

    fit_result_at_E0 = ana.fit_result.parameters.to_table()
    fit_result_at_E0 = fit_result_at_E0[fit_result_at_E0["name"]=="amplitude"]
    fit_result_at_E0["name"] = f"amp@{conf['fit_E0']:5.3}"

    Edec = ana.models[0].spectral_model.pivot_energy.to_value("TeV")
    ana = set_E_ref_to_pivot(ana)
    ana.run_fit()
    print(f"Analytical decorrelation energy is {Edec}")

    fit_result = ana.fit_result.parameters.to_table()
    fit_result = tab.vstack([fit_result,fit_result_at_E0])
    ut.save_fit_result(fit_result, ana, conf, -1.0 * u.TeV)

    best_fit = out_dir / f"{conf['source']}_best_fit.yaml"
    ana.models.write(best_fit, overwrite=True)

    vis.save_fit_residuals(ana.datasets[0], conf, f"PowerLawFit_Residuals")

    if not args.skip_flux:
        #       ana.get_flux_points()
        flux_points, flux_points_est = ut.make_fluxpoints(ana, conf)
        vis.save_flux_points(flux_points_est, conf, f"PowerLawFit_Likelihood")
        try:
            vis.save_fluxpoint_fit(flux_points, conf, f"PowerLawFit")
        except ValueError:
            print("Only upper limits found, residuals can't be calculated")


    opt_res = ana.fit.optimize(ana.datasets)
    ana.fit.covariance(ana.datasets, opt_res)
    fig, cont = vis.plot_contours(conf,ana,"amplitude","index",npoints)
    print(f"Saving countours to: \n"
          f'{conf["out_path"]}/{conf["source"]}_Countours.png')
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_Countours.png')

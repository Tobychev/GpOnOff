#!/usr/bin/env python
import argparse
import os
import pathlib as pt

import gammapy.analysis as ga
import gammapy.data as gda
import gammapy.datasets as gds
import utils as uts
import visualisation as vis
import inputs as inp
import processing as pro

def make_dataset(observations, name, ana_conf, src_pos, conf):
    # ## Setup makers
    (
        ring_bkg_maker,
        safe_mask_maker,
        map_maker,
        excess_maker,
        empty_map,
    ) = pro.setup_no_bkg_ring_makers(ana_conf, src_pos, conf)

    # ## Reduce data

    # full_data = gds.Datasets()
    safe_map = empty_map.copy(name="{name}_stacked")

    print("Doing run:")
    for ob in observations:
        dataset = map_maker.run(empty_map.copy(name=f"{name}_{ob.obs_id}"), ob)
        print(ob.obs_id, end="... ", flush=True)
        try:
            safe_data = safe_mask_maker.run(dataset, ob)
        except:
            print(f"skipping {ob.obs_id}")
            continue

        #full_data.append(dataset_on_off.copy(name=f"{name}_{ob.obs_id}_full"))
        tmp = safe_data.copy(name=f"tmp_{name}_{ob.obs_id}")

        # Does this really make sense?
        safe_data.background = empty_map.exposure
        safe_data.background.data = tmp.exposure.data
        breakpoint()

        safe_ring = ring_bkg_maker.run(safe_data)
        safe_map.stack(safe_ring)

    flux_maps = excess_maker.run(safe_ring)
    print("done")
    return safe_map, flux_maps, excess_maker

# Start of MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to Run ring-background analysis on gammapy formatted data"
    )
    parser.add_argument(
        "config", type=str, help="Config file for the analysis to execute"
    )
    args = parser.parse_args()

    # # Config
    src_ana_conf, src_pos, conf = pro.make_Analysis_config_from_yaml(args.config)
    inp.check_paths(conf)
    rebin = conf["optional"].get("fit_energy_bins", None)
    root_dir = pt.Path(conf["data_directory"]) / conf["optional"].get("production","")
    out_dir = pt.Path(conf["out_path"])

    if "hap-hd" in conf["data_directory"]:
        obs_dict = inp.get_listed_hap_observations(
            conf["optional"]["runlist"], conf["optional"]["cut_conf"]
        )

    if "hap-fr" in conf["data_directory"]:
        obs_dict = inp.get_listed_hd_fr_observations(
            conf["optional"]["runlist"], conf["optional"]["cut_conf"]
        )

    if "fits_export" in conf["data_directory"]:
        obs_dict = inp.get_all_exported(root_dir)


    map_data = gds.Datasets()

    obs_list, prop, names = inp.flatten_obs_dict(obs_dict)

    flux_maps = []
    for lst, name in zip(obs_list, names):
        ring_map, flux_map, exs_maker = make_dataset(lst, name, src_ana_conf, src_pos, conf)
        breakpoint()
        map_data.append(ring_map)
        flux_maps.append(flux_map)


    out_loc = out_dir / "reduced"
    out_loc.mkdir(parents=True, exist_ok=True)

    map_data.write(out_loc / f"{conf['source']}_map_dataset.yaml", overwrite=True)


import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, join, unique
from astropy.time import Time
from gammapy.data import DataStore, get_irfs_features
from gammapy.makers import SafeMaskMaker
from gammapy.maps import WcsGeom
from gammapy.utils.cluster import standard_scaler

HAP_ZEN_BINS = [0, 10, 20, 30, 40, 45, 50, 55, 60, 63, 65, 67, 69, 70]
MUON_LABELS = [
    100,
    101,
    102,
    104,
    105,
    101,
    199,
    200,
    201,
    202,
    203,
    300,
    301,
    302,
    401,
]
MUON_BINS = [
    (20217, 20912),
    (21037, 39875),
    (40093, 56881),
    (57532, 60268),
    (60838, 61615),
    (63567, 67505),
    (68545, 83489),
    (80000, 85000),
    (85000, 95000),
    (95003, 100797),
    (100800, 110340),
    (110340, 124680),
    (127700, 128600),
    (128600, 132350),
    (132350, 154814),
    (154814, 99999999),
]


class DummyDataset:
    def __init__(self, geom):
        self._geom = geom
        self.exposure = None


def fetch_irf_info(run_ids, src_pos, e_ref, store):
    fov_size = 4.0
    thres_meth = "aeff-max"
    thres_meth_par = 10
    runs = store.get_observations(run_ids, required_irf="point-like")
    thresh = np.zeros((len(runs), 2))
    for idx, run in enumerate(runs):
        energy = run.aeff.axes["energy_true"].copy(name="energy")
        geom = WcsGeom.create(
            npix=(1, 1), binsz=fov_size, axes=[energy], skydir=src_pos, proj="TAN"
        )
        dd = DummyDataset(geom)
        smm = SafeMaskMaker(
            methods=[thres_meth, "offset-max"],
            aeff_percent=thres_meth_par,
            position=src_pos,
            irfs="DL3",
        )
        wcsmap = smm.make_mask_energy_aeff_max(dd, run)
        thres_idx = np.where(wcsmap.data)[0][0]

        thresh[idx, 0] = run.obs_id
        thresh[idx, 1] = energy.edges[thres_idx].to("TeV").value

    irf_feats = get_irfs_features(runs, energy_true=e_ref, position=src_pos)
    irf_feats.rename_column("obs_id", "OBS_ID")
    obs_tab = store.obs_table.select_obs_id(run_ids)["OBS_ID", "ZEN_PNT", "ONTIME",
        "LIVETIME", "TSTART", "TSTOP", "DATE-OBS", "TELLIST", "N_TELS", "MUONEFF", "EVENT_COUNT",
        "MUONCORR", "BKG_SCALE"]
    obs_tab = join(obs_tab, Table(
                                thresh,
                             names=["OBS_ID", "E_THR"]),
                    keys="OBS_ID", join_type="left")
    obs_tab = join(obs_tab, irf_feats, keys="OBS_ID", join_type="left")
    obs_tab["MUONEPOC"] = 1
    for lab, (dwn, up) in zip(MUON_LABELS, MUON_BINS):
        sel = (obs_tab["OBS_ID"] >= dwn) & (obs_tab["OBS_ID"] <= up)
        obs_tab["MUONEPOC"][sel] = lab

    return obs_tab


def get_time_groups(obs_tab):
    refi = obs_tab.meta["MJDREFI"]
    reff = obs_tab.meta["MJDREFF"]
    ref_t = Time(refi + reff, format="mjd")

    obs_tab["contigious"] = np.floor(
        obs_tab["TSTOP"] / (50 * u.min.to("s")), dtype=np.int32, casting="unsafe"
    )
    obs_tab["day"] = np.floor(
        obs_tab["TSTOP"] / (1 * u.day.to("s")), dtype=np.int32, casting="unsafe"
    )
    obs_tab["week"] = np.floor(
        obs_tab["TSTOP"] / (1 * u.week.to("s")), dtype=np.int32, casting="unsafe"
    )

    contig = obs_tab.group_by("contigious")
    short_times = []
    for group in contig.groups:
        tsta = group["TSTART"][0] * u.s.to("day") + ref_t
        tsto = group["TSTOP"][-1] * u.s.to("day") + ref_t

        short_times.append(Time([tsta, tsto]))

    contig = obs_tab.group_by("day")
    day_times = []
    for group in contig.groups:
        tsta = group["TSTART"][0] * u.s.to("day") + ref_t
        tsto = group["TSTOP"][-1] * u.s.to("day") + ref_t

        day_times.append(Time([tsta, tsto]))

    contig = obs_tab.group_by("week")
    week_times = []
    for group in contig.groups:
        tsta = group["TSTART"][0] * u.s.to("day") + ref_t
        tsto = group["TSTOP"][-1] * u.s.to("day") + ref_t

        week_times.append(Time([tsta, tsto]))

    return {"contigious": short_times, "daily": day_times, "weekly": week_times}


def make_muon_zenith_plot(obs_table):
    fig, ax = plt.subplots()
    zen = obs_tab["ZEN_PNT"]
    muons = obs_tab["MUONEFF"]
    epoch = obs_tab["MUONEPOC"]
    for epo in np.unique(epoch):
        sel = epoch == epo
        ax.scatter(zen[sel], muons[sel], label=epo)
    ax.legend()
    ax.set(xlabel="Zenith", ylabel="Muon eff")

    return fig


def make_zenith_hist(obs_table):
    zen = obs_tab["ZEN_PNT"]
    fig, ax = plt.subplots()
    ax.hist(zen, bins=HAP_ZEN_BINS)
    ax.set(xlabel="Zenith", ylabel="Frequency")


def make_irf_feat_scatterplots(obs_tab, e_ref="?", axs=None, label=None):
    zen = obs_tab["ZEN_PNT"]
    muons = obs_tab["MUONEFF"]
    ntels = obs_tab["N_TELS"]
    thre = obs_tab["E_THR"]
    bias = obs_tab["edisp-bias"]
    res = obs_tab["edisp-res"]
    psf = obs_tab["psf-radius"]
    rate = obs_tab["EVENT_COUNT"] / obs_tab["LIVETIME"]

    if axs is None:
        fig, axs = plt.subplots(3, 4, figsize=(14, 10), layout="tight")
    elif len(axs.flatten()) < 12:
        raise ValueError("Need at least 12 axes to make scatterplots")

    axs[0, 0].scatter(zen, thre, label=label)
    axs[0, 0].set(xlabel="Zenith", ylabel="Aeff max threshold")

    axs[0, 1].scatter(zen, res, label=label)
    axs[0, 1].set(xlabel="Zenith", ylabel=f"Energy resolution at {e_ref}")

    axs[0, 2].scatter(zen, bias, label=label)
    axs[0, 2].set(xlabel="Zenith", ylabel=f"Energy bias at {e_ref}")

    axs[0, 3].scatter(zen, psf, label=label)
    axs[0, 3].set(xlabel="Zenith", ylabel=f"PSF at {e_ref}")

    axs[1, 0].scatter(muons, thre, label=label)
    axs[1, 0].set(xlabel="Muon eff", ylabel="Aeff max threshold")

    axs[1, 1].scatter(muons, res, label=label)
    axs[1, 1].set(xlabel="Muon eff", ylabel=f"Energy resolution at {e_ref}")

    axs[1, 2].scatter(muons, bias, label=label)
    axs[1, 2].set(xlabel="Muon eff", ylabel=f"Energy bias at {e_ref}")

    axs[1, 3].scatter(muons, psf, label=label)
    axs[1, 3].set(xlabel="Muon eff", ylabel=f"PSF at {e_ref}")

    axs[2, 0].scatter(rate, thre, label=label)
    axs[2, 0].set(xlabel="Event rate", ylabel="Aeff max threshold")

    axs[2, 1].scatter(rate, res, label=label)
    axs[2, 1].set(xlabel="Event rate", ylabel=f"Energy resolution at {e_ref}")

    axs[2, 2].scatter(rate, bias, label=label)
    axs[2, 2].set(xlabel="Event rate", ylabel=f"Energy bias at {e_ref}")

    axs[2, 3].scatter(rate, psf, label=label)
    axs[2, 3].set(xlabel="Event rate", ylabel=f"PSF at {e_ref}")
    return axs

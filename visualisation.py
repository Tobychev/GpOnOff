import gammapy.maps as gmap
import gammapy.visualization as gi
import matplotlib.pyplot as pl


def plot_source_stat(info_table, path=None, prefix='', fig=None, ax_ex=None, ax_sig=None, **kwargs):
    """Plot source statistic for a 1d fit"""
    # make plots on source statistic
    if fig is None:
        fig = pl.figure(figsize = (6,8))
    kwargs.setdefault('ls', 'none')
    kwargs.setdefault('marker', 'o')

    if ax_ex is None or ax_sig is None:
        ax_ex = fig.add_subplot(211)
        ax_sig = fig.add_subplot(212)

    ax_ex.plot(info_table["livetime"].to("h"),
               info_table["excess"],
               **kwargs)

    ax_sig.plot(info_table["livetime"].to("h"),
               info_table["sqrt_ts"],
               **kwargs)

    ax_ex.set_ylabel("Excess")
    ax_sig.set_ylabel("Signficance ($\sigma$)")
    ax_sig.set_xlabel("Lifetime (hours)")
    ax_ex.tick_params(labelbottom = False)

    ax_ex.grid(True)
    ax_sig.grid(True)

    if path is not None:
        fig.subplots_adjust(wspace = 0., hspace = 0., bottom = 0.2, left = 0.2)
        fig.savefig(f"{path}/{prefix}_Source_Stat.png", dpi = 150)
        pl.close("all")

    return ax_ex, ax_sig
    
def plot_spectrum_diagnostic(spectrum):
    fig = pl.figure(figsize=(16,4.5))
    axs = spectrum.peek(fig)
    fig.set_tight_layout(True)
    return fig

def plot_fit_residuals(dataset):
    fig = pl.figure(figsize=(8,6))
    ax1,ax2 = dataset.plot_fit()
    ax1.set_ylabel("Counts")
    ax = dataset.plot_masks(ax1)
    ax.legend()
    return fig

def plot_flux_points(flux_est):
    fig = pl.figure(figsize=(8,6))
    ax = flux_est.plot(sed_type="e2dnde", color="darkorange")
    flux_est.plot_ts_profiles(ax=ax,sed_type="e2dnde")
    return fig

def plot_exp_corrected_counts(spectrum):
   reco_min,reco_max = spectrum.counts.geom.axes["energy"].bounds
   reco_axis = gmap.MapAxis.from_energy_bounds(reco_min.value,reco_max.value,
         nbin = spectrum.counts.geom.axes["energy"].nbin,
         unit = str(reco_min.unit),
         name = "energy_true")
   cnt = spectrum.counts
   exp = spectrum.exposure.resample_axis(reco_axis)
   ene =  cnt.geom.axes.to_table()["ENERGY"].value
   fig = pl.figure(figsize=(8,6))
   pl.loglog(ene,(cnt.data/exp.data).flatten(),'-o')
   return fig

def save_flux_points(flux_est,conf,name):
    fig = plot_flux_points(flux_est)
    fig.suptitle(f"Powerlaw {conf['source']}")
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_{name}.png')
    pl.close("all")

def save_fluxpoint_fit(flux_points,conf,name):
    ax = flux_points.plot_fit()
    fig = pl.gcf()
    fig.suptitle(f"Powerlaw {conf['source']}")
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_{name}.png')
    pl.close("all")

def save_fit_residuals(dataset,conf,name):
    fig = plot_fit_residuals(dataset)

    fig.suptitle(f"Powerlaw {conf['source']} Result")
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_{name}.png')
    pl.close("all")

def save_spectrum_diagnostic(spectrum,conf,name):
    fig = plot_spectrum_diagnostic(spectrum)
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_{name}.png')
    pl.close("all")

def save_off_regions_per_run(dataset,excl_mask):
   for data in datasets:
      pl.figure(figsize(6,6))
      ax = excl_mask.plot()
      print(f"Saving OffRegions{data.name}.png")
      gi.plot_spectrum_datasets_off_regions(ax=ax,datasets=[data])
      pl.savefig(f"OffRegions{data.name}.png")
      pl.close("all")
      
def save_exp_corrected_counts(spectrum,conf,name):
   fig = plot_exp_corrected_counts(spectrum)
   pl.xlabel("Reco energy")
   pl.ylabel("counts/exposure")
   pl.title(conf['source'])
   fig.savefig(f'{conf["out_path"]}/{conf["source"]}_CountsByExp{name}.png')

   pl.close("all")

import gammapy.maps as gmap
import gammapy.visualization as gi
import gammapy.visualization.utils as gvu
import matplotlib.pyplot as pl
import numpy as np

import utils as ut

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
    
def plot_contours(conf,ana,xpar,ypar,npoints=12):
    pars = ut.get_parameter_dict( ana.models.parameters )
    colors = ["m", "b", "c"]
    conts = 3*[None]

    fig,ax = pl.subplots(1,1,figsize=(8,6))
    for sig in (1,2,3):
       idx = sig-1
       try:
          conts[idx] = ana.fit.stat_contour(ana.datasets,
                x=xpar,
                y=ypar,
                sigma=sig,
                numpoints=npoints)
          gvu.plot_contour_line(ax,
               conts[idx][xpar],
               conts[idx][ypar],
               color=colors[idx],
               label=f"{sig}" + r"$\sigma$")
       except:
          print(f"Failed creating countour {sig} using {npoints} points")

    ax.plot(pars[xpar]["value"],
          pars[ypar]["value"],
          marker="x",
          markersize=8)

    if conf["optional"].get("reference_point",False):
       refp = conf["optional"]["reference_point"]

       ax.plot(refp[xpar],
             refp[ypar],
             marker="d",
             markersize=8,
             label=refp["name"])

    ax.set_xlabel(xpar+f" [{pars[xpar]['unit']}]")
    ax.set_ylabel(ypar+f" [{pars[ypar]['unit']}]")
    ax.set_title(f"{conf['source']}")
    pl.legend()
    return fig,conts

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
    # ### Run by run table
    flux_table_path = f'{conf["out_path"]}/{conf["source"]}_{name}_flux.ecsv'
    print(f"Saving flux point table at {flux_table_path}")
    flux_est.to_table().write( flux_table_path ,format="ascii.ecsv", overwrite=True)

def save_fluxpoint_fit(flux_points,conf,name):
    pl.close("all")
    axs = flux_points.plot_fit()
    fig = pl.gcf()
    emin = np.round(
          flux_points.data["energy_min"][0].to("TeV"),
          decimals=1).value
    xlim = axs[0].get_xlim()
    xtics = axs[0].get_xticks()
    sel = xtics > emin

    keep_idx = len(sel) -sum(sel)
    xtics = xtics[(keep_idx-1):]
    xtics[0] =  emin

    xlabs = [str(itm) for itm in xtics]
    xlim = (emin,conf["fit_max"])
    axs[0].set_xticks(xtics)
    axs[1].set_xticks(xtics)
    axs[0].set_xticklabels(xlabs)
    axs[1].set_xticklabels(xlabs)
    axs[0].set_xlim(xlim)
    axs[1].set_xlim(xlim)

    axs[0].text(0.5,0.5,"Preliminary",
          size=50,rotation=30,ha="center",va="center",
          alpha = 0.3,zorder=1000,
          transform=axs[0].transAxes)


    fig.suptitle(f"Powerlaw {conf['source']}")
    fig.savefig(f'{conf["out_path"]}/{conf["source"]}_{name}.png')
    

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


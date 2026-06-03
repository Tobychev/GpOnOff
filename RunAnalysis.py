import hessql
import os
import pandas as pd
import pickle
import numpy as np

from sqlalchemy.exc import OperationalError
pd.set_option('display.max_rows', 500)

def make_runlist_report(runlist,cut_config="spectral",cuts={},short=False,cache=True):

    isfile = False
    if isinstance(runlist,str):
        with open(runlist,"r") as fil:
            runs = [int(itm.split()[0]) for itm in fil.readlines()]
        isfile = True
    elif isinstance(runlist,list):
            runs = runlist
    else:
        raise ValueError(f"runlist must be of type str (a filename), or a list of int. Was {type(runlist)}")

    kind = pd.DataFrame([(itm,get_muon_phase(itm)) for itm in runs],columns=["Run","era"])
    hess1_runs = kind.Run[kind.era < 110]
    hess2_runs = kind.Run[kind.era > 110]

    make_cache = False
    cache_exists = False
    if cache and isfile:

        path,base = os.path.dirname(runlist),os.path.basename(runlist)
        cachefile = f"{path}/.{base}_cache.pkl"
        if not os.path.isfile(cachefile):
            make_cache = True
        else:
            try:
                with open(cachefile,"rb") as fil:
                    qual_data = pickle.load(fil)
                for itm in qual_data:
                    if itm[1] == 1:
                        h1_quality_data = itm[0]
                    if itm[1] == 2:
                        h2_quality_data = itm[0]
                cache_exists = True
            except Exception as Err:
                print("Something wrong with the cache, reloading data from Heidelberg")
                make_cache = True

    quals = []
    if not cache_exists:
        if len(hess1_runs) > 0:
            hess1_db =  hessql.HESS_Database(hessql.hess_database_uri("HD_Monitor"))
        if len(hess2_runs) > 0:
            hess2_db =  hessql.HESS_Database(hessql.hess_database_uri("HD_test"))

        try:
            if len(hess1_runs) > 0:
                h1_quality_data = hess1_db.get_runs_quality(hess1_runs.to_list())
                quals.append((h1_quality_data,1))
            if len(hess2_runs) > 0:
                h2_quality_data = hess2_db.get_runs_quality(hess2_runs.to_list())
                quals.append((h2_quality_data,2))
        except OperationalError as err:
            print("Connection died, please retry")

    if make_cache:
        print(f"Saving quality data to cache file {cachefile}")
        with open(cachefile,"wb") as fil:
            pickle.dump(quals,fil)

    mngrs = []
    if len(hess1_runs) > 0:
        mngrs.append(RunSelectionManager(hess1_runs.to_list(),h1_quality_data,"hess1",cut_config))
    if len(hess2_runs) > 0:
        mngrs.append(RunSelectionManager(hess2_runs.to_list(),h2_quality_data,"hess2",cut_config))

    if len(cuts) > 0:
        for key in cuts.keys():
            for mngr in mngrs:
                mngr.update_cut_criteria(key,cuts[key])

    for mngr in mngrs:
        mngr.make_cut_report()

    if not short:
        for mngr in mngrs:
            if len(mngr.reject_reasons) > 0:
                mngr.details_of_rejected_runs()

    return mngrs

HESS1_SPECTRAL_CUTS = {
    "event_header": lambda data: (data.Participation_frac >= 0.4) & (data.Participation_frac <=1),
    "run_data": lambda data: (data.Duration >= 600) & (data.Duration <=7200),
    "run_atmosphere": lambda data: (data.TransparencyCoefficient_mean >= 0.8) & (data.TransparencyCoefficient_mean <=1.2),
    "run_tracking": lambda data: ( (data.RA_Dev_mean >= -0.01667)  & (data.RA_Dev_mean <= 0.01667) &
                                   (data.Dec_Dev_mean >= -0.01667) & (data.Dec_Dev_mean <= 0.01667) &
                                   (data.Az_Dev_rms >= 0) & (data.Az_Dev_rms <=10) &
                                   (data.Alt_Dev_rms >= 0) & (data.Alt_Dev_rms <=10) ),
    "run_pixel": lambda data: ( (data.Num_Hardware >= 0)  & (data.Num_Hardware <= 120) & 
                                (data.Num_HV_Turned_Off >= 0) & (data.Num_HV_Turned_Off <= 50)  ),
    "run_trigger": lambda data: ( (data.Telescope > 0) & (data.Telescope < 6) &
                    (data.True_Rate_Delta_1 >= -30)  & (data.True_Rate_Delta_1 <= 30) & 
                    (data.True_Rate_Delta_2 >= 0) & (data.True_Rate_Delta_2 <= 10)  )}

HESS2_SPECTRAL_CUTS = {
    "event_header": lambda data: (data.Telescope != 5) & (data.Participation_frac >= 0.04) & (data.Participation_frac <=1),
    "event_header_ct5": lambda data: ((data.Telescope == 5) & 
                                  (data.Participation_frac >= 0.5) & 
                                  (data.Participation_frac <=1) ),
    "event_header_noct5": lambda data: (data.Participation_frac >= 0.4) & (data.Participation_frac <=1) & (data.Telescope != 5),
    "run_data": lambda data: (data.Duration >= 600) & (data.Duration <=7200),
    "run_atmosphere": lambda data: (data.TransparencyCoefficient_mean >= 0.8) & (data.TransparencyCoefficient_mean <=1.2),
    "run_tracking": lambda data: ( (data.RA_Dev_mean >= -0.01667)  & (data.RA_Dev_mean <= 0.01667) & 
                                   (data.Dec_Dev_mean >= -0.01667) & (data.Dec_Dev_mean <= 0.01667) & 
                                   (data.Az_Dev_rms >= 0) & (data.Az_Dev_rms <=10) & 
                                   (data.Alt_Dev_rms >= 0) & (data.Alt_Dev_rms <=10) ),
    "run_pixel": lambda data: ( (data.Num_Hardware >= 0)  & (data.Num_Hardware <= 120) & 
                                (data.Num_HV_Turned_Off >= 0) & (data.Num_HV_Turned_Off <= 50)  ),
    "run_pixel_ct5": lambda data: ( (data.Telescope == 5) &
                                    (data.Num_Hardware >= 0)  & (data.Num_Hardware <= 150) & 
                                (data.Num_HV_Turned_Off >= 0) & (data.Num_HV_Turned_Off <= 80)  ),
    "run_trigger": lambda data: ( (data.Telescope > 0) & (data.Telescope < 6) &
                    (data.True_Rate_Delta_1 >= -30)  & (data.True_Rate_Delta_1 <= 30) & 
                    (data.True_Rate_Delta_2 >= 0) & (data.True_Rate_Delta_2 <= 10)  )}

class RunSelectionManager:

    def __init__(self,runlist,quality_data,era = None, kind = "spectral",mintel = 3, ignore_atmosphere = True, require_CT5 = False):

        self._group_qual = ("event_header","event_header_ct5","event_header_noct5",
                            "run_tracking","run_pixel","run_pixel_ct5","run_trigger")

        if era in ["hess1","hess2"]:
            self.HESS1 = (era == "hess1")
            self.HESS2 = (era == "hess2")
        else:
            raise ValueError("Hess era must have one value out of 'hess1' or 'hess2'")

        if kind == "spectral" and self.HESS1:
            self.cuts = HESS1_SPECTRAL_CUTS
        elif kind == "spectral" and self.HESS2:
            self.cuts = HESS2_SPECTRAL_CUTS
        else:
            raise NotImplemented("f{kind} kind of cuts are not implemted yet")

        self.set_acceptance_critera(mintel, ignore_atmosphere, require_CT5 )

        self.qd = {}
        for key in quality_data:
            self.qd[key] = quality_data[key].copy()


        self.runlist = runlist
        self.reject_reasons = []
        self.unrejected = []

    def details_of_rejected_runs(self,only_no_nans=False):
        if len(self.reject_reasons) < 1:
            raise ValueError("No known reason for rejecting runs, run make_cut_report first")

        if only_no_nans:
            rejected_set = self.value_rejected
        else:
            rejected_set = set(self.missing).union(set(self.value_rejected))

        for reason,count in self.reject_reasons:
            cols = [col for col in self.run_status.columns if reason in col]
            qual_stat = self.run_status.loc[:,cols]
            nas = qual_stat.isna().sum(axis=0).values[0]
            if nas < count:
                qual_stat = qual_stat.loc[(~qual_stat.isna()).sum(axis=1)>0]
            else:
                continue
            print(f"\n\nDetails for criteria {reason}")
            print(f"Out of {count} rejected runs, {count-nas} reject based on {reason} values in database:")
            rejected = ~self.accept[reason](qual_stat)

            tab = self._make_table_name(reason)
            print(self.qd[tab][self.qd[tab].Run.isin(qual_stat.index[rejected])])

    def make_cut_report(self):
        self.run_status = self._place_cuts(self.cuts)
        rejectors = []
        rejects = set()
        nans = 0
        era = "Hess1" if self.HESS1 else "Hess2"
        print(f"\n\nSelection report for {era}:")
        for cond in self.accept.keys():
            sel = ~self.accept[cond](self.run_status)
            if sum(sel) > 0:
                print(f"{sum(sel)} runs discarded by {cond}")
                rejectors.append((cond,sum(sel)))

                rejects = rejects.union(set(sel[sel].index.values))
                if self.run_status[cond].isna().sum() > 0:
                    print(f"   of which {sum(self.run_status[cond].isna())} runs lack info")
                    nans +=1

        if nans > 0:
            nans = self.run_status.isna()
            nansel = nans.sum(axis=1) > 0
            naned = self.run_status.index[nansel]
            nanstat = nans[nansel].sum(axis=0)
            value_rejected = rejects.symmetric_difference(set(naned))
            print(f"\n\nTotal of {len(naned)} runs lack info: ",end="")
            print(set(naned))
            print("Distributed as: ")
            print(nanstat[nanstat > 0])
        else:
            naned = set()
            value_rejected = rejects
        self.missing = naned
        if len(rejects) > 0:
            print(f"\n\nTotal of {len(rejects)} runs discarded ")
            print(rejects,end="\n\n")
            self.reject_reasons = rejectors

        if len(value_rejected) > 0:
            print(f"{len(value_rejected)} runs rejected by value: ",value_rejected )
            self.value_rejected = value_rejected

        else:
            self.value_rejected = set()

        if len(rejects) < len(self.runlist):
            unrejected = rejects.symmetric_difference(set(self.runlist))
            print(f"\n\n{len(unrejected)} runs not rejected: ",unrejected )
            self.unrejected = unrejected

    def update_accept_criteria(self,condition, function):
        self.accept[condition] = function

    def update_cut_criteria(self,condition, function):
        self.cuts[condition] = function

    def set_acceptance_critera(self, mintel = 3, ignore_atmosphere = True, require_CT5 = False):
        self.accept = {
            "participation_quality": lambda df: df.participation_quality >= mintel,
            "duration_quality":lambda df: df.duration_quality,
            "atmosphere_quality":lambda df: df.atmosphere_quality == 6,
            "tracking_quality":lambda df: df.tracking_quality >= mintel,
            "pixel_quality":lambda df: df.pixel_quality >= mintel,
            "trigger_quality":lambda df: df.trigger_quality >= mintel}
        self.mintel = mintel

        if self.HESS2 and require_CT5:
            self.accept["participation_quality"] = lambda df: (
                (df.participation_quality >= (mintel-1))
                 & (df.participation_quality_ct5 == 1))
            self.accept["pixel_quality"] = lambda df: (df.pixel_quality >= (mintel-1)) & (df.pixel_quality_ct5 == 1)
        elif self.HESS2:
            self.accept["participation_quality"] = self._run_by_run_participation_cut
            self.accept["pixel_quality"] = lambda df: (
                ( (df.pixel_quality >= (mintel-1)) & (df.pixel_quality_ct5 == 1)) |
                  (df.pixel_quality >= mintel) )

        if ignore_atmosphere:
             self.accept.pop("atmosphere_quality",None)

    def write_unrejected(self,list_path_name):
        print(f"Writing {len(self.unrejected)} runs to file {list_path_name}")
        with open(list_path_name,"w") as fil:
            fil.writelines([f"{str(runid)}\n" for runid in sorted(self.unrejected)])

    def _run_by_run_participation_cut(self,run_status):
        no_ct5 = lambda df: (df.participation_quality_noct5 >= self.mintel)
        with_ct5 = lambda df: ( (df.participation_quality >= (self.mintel-1)) & (df.participation_quality_ct5 == 1))
        result = run_status["participation_quality"].astype(bool, copy=True)
        rd = self.qd["run_data"]

        for run in run_status.index:
            # This can occasionally fail to identify a 4 tel run
            hess1 = (self.qd["run_data"][ self.qd["run_data"].Run == run].Telescope_Pattern < 31).iloc[0]
            if hess1:
                result[run] = no_ct5(run_status.loc[run])
            else:
                result[run] = with_ct5(run_status.loc[run])
        return result

    def _make_data_name(self,quality_name):
        if "event_header" in quality_name:
            return "event_header"
        elif "run_pixel" in quality_name:
            return "run_pixel"
        else: 
            return quality_name

    def _make_quality_name(self,table_name):
        kind =""
        if "ct5" in table_name:
            kind = "_"+table_name.split("_")[-1]
        if "event_header" in table_name:
            return "participation_quality"+kind
        elif table_name == "run_data":
            return "duration_quality"
        elif "run_" in table_name:
            return "{}_quality".format(table_name.split("_")[1])+kind
        else: 
            raise NotImplementedError

    def _make_table_name(self,quality):
        names = {"participation_quality":"event_header",
                 "duration_quality":"run_data",}
        if quality in names.keys():
            return names[quality]
        else:
            return f"run_{quality.split('_')[0]}"

    def _place_cuts(self,cuts):
        selection = {}
        for key in cuts.keys():
            tab_name = self._make_quality_name(key)
            dat_name = self._make_data_name(key)
            if key in self._group_qual:
                # Evaluate per-telescope cut and insert an indicator for each telescope
                # found. Then aggregate to find the per-run total quality
                if not tab_name in self.qd[dat_name].columns:
                    self.qd[dat_name].insert(2,tab_name,cuts[key](self.qd[dat_name]))
                grouped = self.qd[dat_name].groupby("Run")[tab_name]
                selection[key] = grouped.agg(sum)
                # In additions, if there were not enough telescopes present to 
                # satisfy the condition, mark this with NaN so it is not counted
                # like the run was cut because it was out of range
                avail = grouped.count() >= self.mintel
                selection[key][np.invert(avail)] = np.nan
            elif key == "run_data":
                if not tab_name in self.qd[dat_name].columns:
                    self.qd[dat_name].insert(2,tab_name,cuts[key](self.qd[dat_name]))
                qual_tab = self.qd[dat_name][["Run",tab_name]].copy()
                qual_tab.set_index("Run",inplace=True)
                # Set rows that failed due to missing telescopes to nan
                avail = self.qd[dat_name]["Number_Of_Telescopes"] < self.mintel
                avail = avail.to_frame()
                qual_tab[avail.fillna(True)]
                selection[key] = qual_tab

            elif key == "run_atmosphere":
                if not tab_name in self.qd[dat_name].columns:
                    self.qd[dat_name].insert(2,tab_name,cuts[key](self.qd[dat_name]))
                selection[key] = self.qd[dat_name][["Run",tab_name]].copy()
                selection[key].set_index("Run",inplace=True)
                selection[key] =  selection[dat_name]*6

            else:
                raise ValueError("Uknonwn cut name {}".format(key))

        self.selection = pd.concat([itm for k,itm in selection.items()],axis=1,join="outer")
        return self.selection.copy()


def get_muon_phase(run):

    phase1 = 20000 # Actually the start of reliable runs
    phase1b = 20987
    phase1c = 39998
    phase1c1 = 57434
    phase1c2 = 60742
    phase1c3 = 63295
    phase1d = 68311
    phase2b0 = 79975
    phase2b2 = 84921
    phase2b3 = 95010
    phase2b4 = 100796
    phase2b5 = 110322
    phase2c0 = 127583
    phase2c1 = 128592
    phase2c2 = 131717
    phase3c0 = 154814

    if (run >= phase1 )  & (run < phase1b) : return 100
    if (run >= phase1b)  & (run < phase1c) : return 101
    if (run >= phase1c)  & (run < phase1c1): return 102
    if (run >= phase1c1) & (run < phase1c2): return 103
    if (run >= phase1c2) & (run < phase1c3): return 104
    if (run >= phase1c3) & (run < phase1d) : return 105
    if (run >= phase1d ) & (run < phase2b0): return 101
    if (run >= phase2b0) & (run < phase2b2): return 199
    if (run >= phase2b2) & (run < phase2b3): return 200
    if (run >= phase2b3) & (run < phase2b4): return 201
    if (run >= phase2b4) & (run < phase2b5): return 202
    if (run >= phase2b5) & (run < phase2c0): return 203
    if (run >= phase2c0) & (run < phase2c1): return 300
    if (run >= phase2c1) & (run < phase2c2): return 301
    if (run >= phase2c2) & (run < phase3c0): return 302
    if (run >= phase3c0)                   : return 402        

def short_run_quality_report(runlist,era,kind):
    if not isinstance(runlist,list):
        raise ValueError(f"runlist must be of type list, was {type(runlist)}")
    if era == "hess1":
        hess =  hessql.HESS_Database(hessql.hess_database_uri("HD_Monitor"))
    elif era == "hess2":
        hess =  hessql.HESS_Database(hessql.hess_database_uri("HD_test"))
    try:
        quality_data = hess.get_runs_quality(runlist)
        runs = RunSelectionManager(quality_data,era=era,kind=kind)
    except OperationalError as err:
        print("Connection died, please retry")

    rejectors = []
    rejects = set()
    nans = 0
    print("\n\nSelection report:")
    for cond in runs.accept.keys():
        sel = ~runs.accept[cond](runs.run_status)
        if sum(sel) > 0:
            print(f"{sum(sel)} runs discarded by {cond}")
            rejectors.append((cond,sum(sel)))

            rejects = rejects.union(set(sel[sel].index.values))
            if runs.run_status[cond].isna().sum() > 0:
                print(f"   of which {sum(runs.run_status[cond].isna())} runs lack info")
                nans +=1

    if nans > 0:
        nans = runs.run_status.isna()
        nansel = nans.sum(axis=1) > 0
        naned = runs.run_status.index[nansel]
        nanstat = nans[nansel].sum(axis=0)
        value_rejected = rejects.symmetric_difference(set(naned))
        print(f"\n\nTotal of {len(naned)} runs lack info: ",end="")
        print(set(naned))
        print("Distributed as: ")
        print(nanstat[nanstat > 0])
    else:
        naned = None
        value_rejected = rejects
    if len(rejects) > 0:
        print(f"\n\nTotal of {len(rejects)} runs discarded ")
        print(rejects,end="\n\n")

    if len(value_rejected) > 0:
        print(f"{len(value_rejected)} runs rejected by value: ",value_rejected )

    unrejected = rejects.symmetric_difference(set(runlist))
    if len(unrejected) > 0:
        print(f"{len(unrejected)} runs not rejected: ",unrejected )

    return {"rejected":rejects,"selector":runs,"reasons":rejectors,"missing":naned}

def details_of_rejected_runs(ra_result,only_no_nans=False):
    tot_rejected = ra_result["rejected"]
    missing = ra_result["missing"]
    sels = ra_result["selector"]
    qual_rejected = missing.symmetric_difference(set(tot_rejected))
    if only_no_nans:
        rejected_set = qual_rejected
    else:
        rejected_set = tot_rejected

    for reason,count in ra_result["reasons"]:
        cols = [col for col in sels.run_status.columns if reason in col]
        qual_stat = sels.run_status.loc[:,cols]
        nas = qual_stat.isna().sum(axis=0).values[0]
        print(nas)
        if nas < count:
            qual_stat = qual_stat.loc[(~qual_stat.isna()).sum(axis=1)>0]
        else:
            continue
        print(f"\n\nOut of {count} rejected runs, {count-nas} reject based on {reason} values in database:")
        rejected = ~sels.accept[reason](qual_stat)

        tab = sels._make_table_name(reason)
        print(sels.qd[tab][sels.qd[tab].Run.isin(qual_stat.index[rejected])])

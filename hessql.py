# +
import sqlalchemy as sqla
import configparser
import copy
import pandas as pd

from sqlalchemy.sql.expression import func
# -

import astropy.coordinates as ac

from collections import namedtuple as nt

database = {"hess1":"HD_Monitor","hess2":"HD_Test"}

def hess_database_uri(base): 
    conf = configparser.ConfigParser()
    conf.read("/nfs/us0/tobychev/.dbtoolsrc")
    return 'mysql+pymysql://{user}:{password}@{host}:{port}/{database}?charset=utf8'.format(
    host     = conf["hess"]['host'],
    user     = conf["hess"]['user'],
    password = conf["hess"]['password'],
    port     = conf["hess"]['port'],
    database = base)

class HESS_Database:

    def __init__(self,database_uri):
        self.DATABASE_URI = database_uri
        self.engine = sqla.create_engine(self.DATABASE_URI)
        self.meta = sqla.MetaData()
        self.meta.reflect(bind=self.engine)

        self.evn_head = self.meta.tables['Monitor_Run_EventHeader']
        self.run_data = self.meta.tables['Monitor_Run_Data']
        self.run_trak = self.meta.tables['Monitor_Run_Tracking']
        self.run_pixl = self.meta.tables['Monitor_Run_Pixel']
        self.run_trig = self.meta.tables['Monitor_Run_Trigger']
        self.run_atmo = self.meta.tables['Monitor_Run_Atmosphere']

        self.qualities = {
            "event_header":self.evn_head,
            "run_data":self.run_data,
            "run_tracking":self.run_trak,
            "run_pixel":self.run_pixl,
            "run_trigger":self.run_trig,
            "run_atmosphere":self.run_atmo
        }

    def make_latest_or(self,table,runlist):
        latest_rows = sqla.select(
            table.c["Run"],
            func.max(table.c["WhenEntered"])
        ).where(
            table.c["Run"].in_(runlist)
            ).group_by(table.c["Run"])

        latest_qual = pd.read_sql_query(latest_rows,self.engine)
        # Some rows have been updated twice the same day, lets
        # hope they haven't been updated twice within 30 seconds
        return [sqla.and_(table.c["Run"] == run,sqla.between(table.c["WhenEntered"],
             (date-pd.Timedelta(seconds=30)).to_pydatetime(), (date+pd.Timedelta(seconds=2)).to_pydatetime() 
                                                     )
                  )
         for run,date in latest_qual.values]

    def select_event_quality(self,latest_or_condition):
        hed_qual = sqla.select(
            self.evn_head.c["Run"],
            self.evn_head.c["Telescope"],
            self.evn_head.c["Participation_frac"]
        ).where(
            sqla.and_(
                self.evn_head.c["Telescope"] != 0,
                sqla.or_(*latest_or_condition) )
            )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_data_quality(self,latest_or_condition):
        hed_qual = sqla.select(
            self.run_data.c["Run"],
            self.run_data.c["Duration"],
            self.run_data.c["Offset_x"],
            self.run_data.c["Offset_y"],
            self.run_data.c["Telescope_Pattern"],
            self.run_data.c["Number_Of_Telescopes"],
            self.run_data.c["RunType"],
            self.run_data.c["Offset_Mode"]
        ).where(
                sqla.or_(*latest_or_condition) )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_tracking_quality(self,latest_or_condition):
        hed_qual = sqla.select(
            self.run_trak.c["Run"],
            self.run_trak.c["Telescope"],
            self.run_trak.c["RA_Dev_mean"],
            self.run_trak.c["Dec_Dev_mean"],
            self.run_trak.c["Az_Dev_rms"],
            self.run_trak.c["Alt_Dev_rms"], 
        ).where(
            sqla.and_(
                self.run_trak.c["Telescope"] != 0,
                sqla.or_(*latest_or_condition) )
            )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_pixel_quality(self,latest_or_condition):
        hed_qual = sqla.select(
            self.run_pixl.c["Run"],
            self.run_pixl.c["Telescope"],
            self.run_pixl.c["Num_Hardware"],
            self.run_pixl.c["Num_Broken"],
            self.run_pixl.c["Num_HV_Turned_Off"],
        ).where(
            sqla.and_(
                self.run_pixl.c["Telescope"] != 0,
                sqla.or_(*latest_or_condition) )
            )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_atmosphere_quality(self,latest_or_condition):
        hed_qual = sqla.select(
            self.run_atmo.c["Run"],
            self.run_atmo.c["TransparencyCoefficient_mean"],
        ).where(
                sqla.or_(*latest_or_condition) )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_tigger_quality(self,latest_or_condition):
        columns = [ self.run_trig.c["Run"],
                 self.run_trig.c["Telescope"],
                 self.run_trig.c["Mean_Zenith"],
                 self.run_trig.c["True_Rate_mean"],
                 self.run_trig.c["True_Rate_Delta_1"],
                 self.run_trig.c["True_Rate_Delta_2"]]

        hed_qual = sqla.select( *columns
        ).where(
            sqla.and_(
                self.run_trig.c["Telescope"] != 0,
                self.run_trig.c["Telescope"] < 6,
                sqla.or_(*latest_or_condition) )
            )
        return pd.read_sql_query(hed_qual,self.engine)

    def select_quality(self,quality,or_condition):
        if quality == "event_header":
            return self.select_event_quality(or_condition)
        if quality == "run_data":
            return self.select_data_quality(or_condition)
        if quality == "run_tracking":
            return self.select_tracking_quality(or_condition)
        if quality == "run_pixel":
            return self.select_pixel_quality(or_condition)
        if quality == "run_trigger":
            return self.select_tigger_quality(or_condition)
        if quality == "run_atmosphere":
            return self.select_atmosphere_quality(or_condition)
        else:
            return None

    # to get the raw sql: print(latest_rows.compile(compile_kwargs={"literal_binds": True}))
    def get_runs_quality(self,runlist):
        data = {}
        for qual in self.qualities:
            print("Fetching quality data for {}".format(qual))
            or_condition = self.make_latest_or(self.qualities[qual],runlist)
            if len(or_condition) > 0:
                data[qual] = self.select_quality(qual,or_condition)
            else:
                # If there aren't any latest runs ensure we return nothing
                # or_condition has to be iterable
                data[qual] = self.select_quality(qual,(False,False))

        return data

    def get_runs_in_coord_box(self, coords, margin, only_obs=True):
        if only_obs:
            runs = sqla.select(self.run_data).\
                                where(
                                    sqla.and_(
                                        self.run_data.c['RunType'] == "ObservationRun" ,
                                        sqla.and_(
                                            self.run_data.c["Target_RA"] > coords.ra.value - margin,
                                            self.run_data.c["Target_RA"] < coords.ra.value + margin),
                                        sqla.and_(
                                            self.run_data.c["Target_Dec"] > coords.dec.value - margin,
                                            self.run_data.c["Target_Dec"] < coords.dec.value + margin),
                                    )
                                )
        else:
            runs = sqla.select(self.run_data).\
                                where(
                                    sqla.and_(
                                        sqla.and_(
                                            self.run_data.c["Target_RA"] > coords.ra.value - margin,
                                            self.run_data.c["Target_RA"] < coords.ra.value + margin),
                                        sqla.and_(
                                            self.run_data.c["Target_Dec"] > coords.dec.value - margin,
                                            self.run_data.c["Target_Dec"] < coords.dec.value + margin),
                                    )
                                )

        return pd.read_sql_query(runs,self.engine)

class SelectQualityRuns:

    def __init__(self,quality_data):
        self._group_qual = ("event_header","event_header_ct5","run_tracking","run_pixel","run_pixel_ct5","run_trigger")

        self.qd = {}
        for key in quality_data:
            self.qd[key] = quality_data[key].copy()

    def make_data_name(self,quality_name):
        if "event_header" in quality_name:
            return "event_header"
        elif "run_pixel" in quality_name:
            return "run_pixel"
        else: 
            return quality_name

    def make_quality_name(self,table_name):
        if "event_header" in table_name:
            return "participation_quality"
        elif table_name == "run_data":
            return "duration_quality"
        elif "run_" in table_name:
            return "{}_quality".format(table_name.split("_")[1])
        else: 
            raise NotImplementedError

    def place_cuts(self,cuts):
        selection = {}
        for key in cuts.keys():
            tab_name = self.make_quality_name(key)
            if "ct5" in key:
                tab_name += "_ct5"
            dat_name = self.make_data_name(key)
            if key in self._group_qual:
                if not tab_name in self.qd[dat_name].columns:
                    self.qd[dat_name].insert(2,tab_name,cuts[key](self.qd[dat_name]))
                selection[key] = self.qd[dat_name].groupby("Run")[tab_name].agg(sum)

            elif key == "run_data":
                if not tab_name in self.qd[dat_name].columns:
                    self.qd[dat_name].insert(2,tab_name,cuts[key](self.qd[dat_name]))
                selection[key] = self.qd[dat_name][["Run","Number_Of_Telescopes"]].copy()
                selection[key].set_index("Run",inplace=True)
                selection[key].rename(inplace=True,
                    columns={"Number_Of_Telescopes":tab_name})

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


    def select_passing(self,num_of_tel = 3, ct5 = None, ignore = None ,selection = None):
        if not selection:
            sel = self.selection
        else:
            sel = selection

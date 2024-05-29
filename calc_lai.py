# -*- coding: utf-8 -*-
# Copyright (c) 2024 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), May 2024
# ===================================================
"""
Sample script loading a parcel and computing LAI based on weather data
downloaded with download_weather.py and the LAI generator class in
utils/lai_generator.py. More info in both scripts.
"""
import warnings

from agromanagement.utility.lai_generator import LaiGenerator
from agromanagement.utility.paths import ROOT_DIR

warnings.simplefilter(action="ignore")

PARCEL_NAME = "723134"
parcel_path = f"D:/Documents/GitHub/AgroManagement/resources/{PARCEL_NAME}.geojson"

SENTINEL_1_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/"
SENTINEL_2_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/"
ERA_DIR = "resources/era_5/"

START_DATE = "2019-01-01"
END_DATE = "2019-12-31"

lai = LaiGenerator(
    jsonloc=parcel_path,
    filename="prova",
    workingdir=ROOT_DIR,
    s1_dir=SENTINEL_1_UNZIPPED_FOLDER,
    s2_dir=SENTINEL_2_UNZIPPED_FOLDER,
    era_dir=ERA_DIR,
    snap_graphs_dir="C:/Users/mcm216/.snap/graphs/",
    snap_gtp_dir="C:/Program Files/esa-snap/bin/",
    startdate=START_DATE,
    enddate=END_DATE,
)

lai.s1_to_vvvh()
lai.s2_to_lai()
lai.lai_ts_creation()

# -*- coding: utf-8 -*-
# Copyright (c) 2024 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), May 2024
# ===================================================
"""
This is a sample script to download the weather data required to
test the script that estimates LAI based on Sentinel 1, Sentinel 2
and ERA 5 reanalysis data.
The script is largely based on the example notebooks 2 and 3 in the
'examples' folder. The main difference is that here we use the
cdsetool.download:download_features() function rather than the
cdsetool.download:download_feature(), which bulk downloads the data
rather than downloading only one feature at the time.
"""
import sys

import geopandas as gpd
from cdsetool.credentials import Credentials
from cdsetool.download import download_features
from cdsetool.query import query_features
from shapely.geometry import mapping

from agromanagement.utility.data_downloaders import download_era
from agromanagement.utility.file_management import unzip_all_files
from agromanagement.utility.lai_fusion import LAI_S1S2_fusion
from agromanagement.utility.paths import ROOT_DIR

sys.path.append('../')

credentials = Credentials()


# define a parcel
parcels = gpd.read_file("resources/" + "lcm2021_tile_11_1014_647.geojson")
parcel = parcels.iloc[1025:1026,:]
print(parcel)

geometry = parcel.iloc[0,:]["geometry"]
geom_json = mapping(geometry)

START_DATE = "2019-04-01"
END_DATE = "2019-04-30"

SENTINEL_1_OUTPUT_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/Raw/"
SENTINEL_2_OUTPUT_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/Raw/"


S1_COLLECTION = "Sentinel1"
S1_PRODUCT_TYPE = "IW_GRDH_1S"

search_terms = {
    "maxRecords": "2000",
    "startDate": START_DATE,
    "completionDate": END_DATE,
    "geometry": geometry,
    "productType": S1_PRODUCT_TYPE,
}

sentinel_1_feature_list = query_features(collection=S1_COLLECTION, search_terms=search_terms)


S2_COLLECTION = "Sentinel2"
S2_PRODUCT_TYPE = "S2MSI2A"

search_terms = {
    "maxRecords": "2000",
    "startDate": START_DATE,
    "completionDate": END_DATE,
    "geometry": geometry,
    "productType": S2_PRODUCT_TYPE,
}

sentinel_2_feature_list = query_features(collection=S2_COLLECTION, search_terms=search_terms)

s1_downloads = download_features(sentinel_1_feature_list, SENTINEL_1_OUTPUT_FOLDER)
for item in s1_downloads:
    print(f"feature {item} downloaded")

s2_downloads = download_features(sentinel_2_feature_list, SENTINEL_2_OUTPUT_FOLDER)
for item in s2_downloads:
    print(f"feature {item} downloaded")

download_era(start_date=START_DATE, end_date=END_DATE)

SENTINEL_1_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/Unzipped/"
SENTINEL_2_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/Unzipped/"

unzip_all_files(input_folder=SENTINEL_1_UNZIPPED_FOLDER)
unzip_all_files(
    input_folder=SENTINEL_2_OUTPUT_FOLDER,
    output_folder=SENTINEL_2_UNZIPPED_FOLDER
)

## LAI attempt
PARCEL_PATH = "resources/" + "lcm2021_tile_11_1014_647.geojson"
ERA_DIR = "resources/era_5/"
lai = LAI_S1S2_fusion(
    jsonloc=PARCEL_PATH,
    FNAME="prova",
    workingdir=ROOT_DIR,
    s1dir=SENTINEL_1_UNZIPPED_FOLDER,
    s2dir=SENTINEL_2_UNZIPPED_FOLDER,
    ECMWF_ERA5_VPD=ERA_DIR,
    snap_graphs_dir="resources/snap_graphs/",
    snap_gtp_dir="C:/Program Files/esa-snap/bin/",
    startdate=START_DATE,
    enddate=END_DATE
)

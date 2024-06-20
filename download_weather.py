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
from cdsetool.monitor import StatusMonitor
from cdsetool.query import query_features
from shapely.geometry import mapping

from agromanagement.utility.data_downloaders import download_era
from agromanagement.utility.file_management import unzip_all_files

sys.path.append("../")

credentials = Credentials()

# define a parcel
parcel = gpd.read_file("resources/" + "SY17219007.geojson")

geometry = parcel.iloc[0, :]["geometry"]
geom_json = mapping(geometry)

START_DATE = "2019-01-01"
END_DATE = "2019-12-31"

SENTINEL_1_OUTPUT_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/SY17219007/"
SENTINEL_2_OUTPUT_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/SY17219007/"


S1_COLLECTION = "Sentinel1"
S1_PRODUCT_TYPE = "IW_GRDH_1S"

search_terms = {
    "maxRecords": "2000",
    "startDate": START_DATE,
    "completionDate": END_DATE,
    "geometry": geometry,
    "productType": S1_PRODUCT_TYPE,
}

sentinel_1_feature_list = query_features(
    collection=S1_COLLECTION, search_terms=search_terms
)

S2_COLLECTION = "Sentinel2"
S2_PRODUCT_TYPE = "S2MSI2A"

search_terms = {
    "maxRecords": "2000",
    "startDate": START_DATE,
    "completionDate": END_DATE,
    "geometry": geometry,
    "productType": S2_PRODUCT_TYPE,
    "processingBaseline": "05.00",
}

sentinel_2_feature_list = query_features(
    collection=S2_COLLECTION, search_terms=search_terms
)
len(sentinel_2_feature_list)

list(
    download_features(
        sentinel_1_feature_list,
        SENTINEL_1_OUTPUT_FOLDER,
        {
            "concurrency": 5,
            "monitor": StatusMonitor(),
            "credentials": Credentials(),
        },
    )
)

list(download_features(
        sentinel_2_feature_list,
        SENTINEL_2_OUTPUT_FOLDER,
        {
            "concurrency": 5,
            "monitor": StatusMonitor(),
            "credentials": Credentials(),
        },
    )
)

download_era(start_date=START_DATE, end_date=END_DATE)

## Unizip if required
SENTINEL_1_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/Unzipped/"
SENTINEL_2_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/Unzipped/"

unzip_all_files(input_folder=SENTINEL_1_UNZIPPED_FOLDER)
unzip_all_files(
    input_folder=SENTINEL_2_OUTPUT_FOLDER, output_folder=SENTINEL_2_UNZIPPED_FOLDER
)

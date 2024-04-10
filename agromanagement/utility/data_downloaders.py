# -*- coding: utf-8 -*-
# Copyright (c) 2023 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), September 2023
# =========================================================
"""
Utility functions to download and store a variety of climate data.
==================================================================

Functions defined here:
-----------------------

download_sentinel(geotiff_file, output_folder):
    Download a geotiff file from Sentinel 2 L2A after cropping to match 
    the bounding box defined in bbox.

"""
import os
from datetime import datetime

import cdsapi
import xarray as xr

# import rasterio
# from agromanagement.core import aws_session
from agromanagement.utility.paths import ROOT_DIR

# def download_sentinel(geotiff_file, output_folder=None):
#     """
#     Download a geotiff file from Sentinel 2 L2A after cropping to match
#     the bounding box defined in bbox. This works only for single band tiff
#     files. The geotiff file is saved in output_folder/tile_subfolder
#     where the 'tile_subfolder' is created based on the name of the tile
#     of interest. If looping through various bands, this function will save
#     all bands in the same subfolder.

#     Parameters
#     ----------
#     :param geotiff_file (str): the link ['href'] of the location of the
#         geotiff file of interest. For example:
#         >>> sentinel_items = SentinelSearch.items()
#         >>> item = sentinel_items[0]
#         >>> geotiff_file = item.assets["B04"]["href2]
#         (-2.116652, 60.139989, -2.116526, 60.14011)
#     :param output_folder (str), optional: the path where the data will be saved.
#         If not none, the function will create inside output_folder a subfolder
#         (if missing) whose name is retrieved from 'geotiff_file' based on
#         standard Sentinel 2 naming conventions, and the data will be saved
#         inside that subfolder. This makes sure that if multiple 'geotiff_file'
#         (i.e., multiple bands for the same tile) are downloaded, for example
#         in a loop calling this function, all bands from the same tile are
#         saved in the same subfolder.
#     : param save_file (bool): whether to save th
#     Return
#     ------

#     """
#     with rasterio.Env(aws_session):
#         with rasterio.open(geotiff_file) as geo_fp:
#             folder_name = geotiff_file[-32:-8]
#             band_name = geotiff_file[-7:]
#             if output_folder is None:
#                 metadata = geo_fp.meta
#                 metadata.update(
#                     {
#                         "driver": "GTiff",
#                         "count": 1,  # Assuming single band
#                     }
#                 )
#                 tile_data = geo_fp.read(1)
#                 print(
#                     f"File '{band_name}' for tile {folder_name}"
#                     f" successfully downloaded."
#                 )
#                 return {"data": tile_data, "metadata": metadata}
#             tile_folder = os.path.join(output_folder, folder_name)
#             if not os.path.exists(tile_folder):
#                 os.makedirs(tile_folder)
#             output_filename = f"{tile_folder}\\{band_name}"
#             if os.path.exists(output_filename):
#                 raise FileExistsError(
#                     f"The file '{band_name}' for tile {folder_name}" f" already exists."
#                 )
#             meta = geo_fp.meta
#             meta.update(
#                 {
#                     "driver": "GTiff",
#                     "count": 1,  # Assuming single band
#                 }
#             )
#             tile_data = geo_fp.read(1)
#             with rasterio.open(output_filename, "w", **meta) as dst:
#                 dst.write(tile_data, 1)
#             print(
#                 f"File '{band_name}' for tile {folder_name}"
#                 f" successfully downloaded."
#             )
#             return None

DEFAULT_DOWNLOAD_PATH = os.path.join(ROOT_DIR, "resources", "era_5")
if not os.path.exists(DEFAULT_DOWNLOAD_PATH):
    os.makedirs(DEFAULT_DOWNLOAD_PATH)


def download_era(start_date, end_date, download_path=DEFAULT_DOWNLOAD_PATH):
    """
    Download and store a netcdf file containing ERA5 reanalysis
    data to compute Vapour Pressure Deficit to use to estimate
    LAI in the UK.
    We use Copernicus and we interface with their servers using
    the CDS API:
    (https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5).
    More info on how to optimise data download can be found here:
    http://tinyurl.com/5dvy4evm

    Parameters
    ----------
    :param start_date (str): The start date for the timeframe of
        interest. Declared as "yyyy-mm-dd", e.g., "2019-01-01"
    :param end_date (str): The end date for the timeframe of
        interest. Declared as "yyyy-mm-dd", e.g., "2021-12-31"
    :param download_path (str): the location where the downloaded
        data will be stored. Default set to the /resources/era_5
        folder in the main project directory.

    N.B.: the function will always retrieve and save reanalysis data
        for the entire years in the range between start_date and
        end_date, regardless of whether the timeframe only covers part
        of the years.
    """
    cds_client = cdsapi.Client()
    start_datetime = datetime.strptime(start_date, "%Y-%m-%d").date()
    end_datetime = datetime.strptime(end_date, "%Y-%m-%d").date()
    for year in range(start_datetime.year, end_datetime.year + 1):
        yearly_filename = f"{download_path}/ERA5_{year}.nc"
        if os.path.exists(yearly_filename):
            print(f"Data for year '{year}' already downloaded. Skipping...")
            continue
        file_list = []
        for month in range(1, 13):
            print("========================================================")
            print(f"Downloading data for year '{year}' and month '{month}' ...")
            monthly_filename = f"{download_path}/ERA5_{year}_{month:02d}.nc"
            if os.path.exists(monthly_filename):
                print(
                    f"data for year '{year}' and month '{month}' "
                    f"already exists. Skipping..."
                )
                file_list.append(monthly_filename)
                continue
            cds_client.retrieve(
                "reanalysis-era5-land",
                {
                    "product_type": "reanalysis",
                    "format": "netcdf",
                    "variable": [
                        "2m_temperature",
                        "2m_dewpoint_temperature",
                    ],
                    "year": str(year),
                    "month": f"{month:02d}",
                    "day": [str(i).zfill(2) for i in range(1, 32)],
                    "time": [f"{hour:02d}:00" for hour in range(24)],
                    "area": [61, -9, 49, 2],
                },
                f"{download_path}/ERA5_{year}_{month:02d}.nc",
            )
            file_list.append(monthly_filename)

        if len(file_list) != 12:
            raise ValueError(
                f"There are not 12 months of available data for year '{year}'."
            )
        print(f"Combining monthly data for year '{year} ...'")
        datasets = [xr.open_dataset(file) for file in file_list]
        combined_dataset = xr.concat(datasets, dim="time")
        combined_dataset.to_netcdf(f"{download_path}/ERA5_{year}.nc")

        for dataset in datasets:
            dataset.close()
            del dataset

        for file in file_list:
            os.remove(file)
        print("... done ...")

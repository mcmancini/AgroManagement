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

# import os
# import rasterio
# from agromanagement.core import aws_session
# from agromanagement.utility.paths import ROOT_DIR

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

# DEFAULT_DOWNLOAD_PATH = os.path.join(ROOT_DIR, "resources", "era_5")
# def download_era(start_date, end_date, download_path):
#     """
#     Download and store a netcdf file containing ERA5 reanalysis
#     data to compute Vapour Pressure Deficit to use to estimate
#     LAI in the UK.
#     We use Copernicus and we interface with their servers using
#     the CDS API:
#     (https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5).
#     More info on how to optimise data download can be found here:
#     http://tinyurl.com/5dvy4evm

#     Parameters
#     ----------
#     :param start_date (str): The start date for the timeframe of
#         interest. Declared as "yyyy-mm-dd", e.g., "2019-01-01"
#     :param end_date (str): The end date for the timeframe of
#         interest. Declared as "yyyy-mm-dd", e.g., "2021-12-31"
#     :param
#     """
#     pass

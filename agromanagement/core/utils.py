# -*- coding: utf-8 -*-
# Copyright (c) 2023 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), September 2023
# =========================================================
"""
Utility functions to be used in the UK implementation of the
WOFOST crop yield model.
============================================================

Functions defined here:
-----------------------

lonlat2osgrid(coords, figs)
    converts a lon-lat pair to an OS grid code

osgrid2lonlat(gridref, epsg)
    Convert British National Grid references to OSGB36 numeric
    coordinates. Grid references can be 4, 6, 8 or 10 figures

nearest(item, valuelist):
    Find nearest value to item in valuelist

get_dtm_values(parcel_os_code, app_config)
    Query the DTM database based on longitude and latitude to
    retrieve elevation, slope and aspect

read_parcel_data(gid, folder):
    Read UKCEH parcel data retrieved from a series of
    files contained in "folder".

download_sentinel(geotiff_file, output_folder):
    Download a geotiff file from Sentinel 2 L2A after cropping to match 
    the bounding box defined in bbox.
"""
import json
import os
import re
from math import floor

import psycopg2
from pyproj import Transformer


class BNGError(Exception):
    """Exception raised by OSgrid coordinate conversion functions"""


def _init_regions_and_offsets():
    # Region codes for 100 km grid squares.
    regions = [
        ["HL", "HM", "HN", "HO", "HP", "JL", "JM"],
        ["HQ", "HR", "HS", "HT", "HU", "JQ", "JR"],
        ["HV", "HW", "HX", "HY", "HZ", "JV", "JW"],
        ["NA", "NB", "NC", "ND", "NE", "OA", "OB"],
        ["NF", "NG", "NH", "NJ", "NK", "OF", "OG"],
        ["NL", "NM", "NN", "NO", "NP", "OL", "OM"],
        ["NQ", "NR", "NS", "NT", "NU", "OQ", "OR"],
        ["NV", "NW", "NX", "NY", "NZ", "OV", "OW"],
        ["SA", "SB", "SC", "SD", "SE", "TA", "TB"],
        ["SF", "SG", "SH", "SJ", "SK", "TF", "TG"],
        ["SL", "SM", "SN", "SO", "SP", "TL", "TM"],
        ["SQ", "SR", "SS", "ST", "SU", "TQ", "TR"],
        ["SV", "SW", "SX", "SY", "SZ", "TV", "TW"],
    ]

    # Transpose so that index corresponds to offset
    regions = list(zip(*regions[::-1]))

    # Create mapping to access offsets from region codes
    offset_map = {}
    for i, row in enumerate(regions):
        for j, region in enumerate(row):
            offset_map[region] = (1e5 * i, 1e5 * j)

    return regions, offset_map


_regions, _offset_map = _init_regions_and_offsets()


def lonlat2osgrid(coords, figs=4):
    """
    Convert WGS84 lon-lat coordinates to British National Grid references.
    Grid references can be 4, 6, 8 or 10 fig, specified by the figs keyword.
    Adapted from John A. Stevenson's 'bng' package that can be found at
    https://pypi.org/project/bng/

    :param coords: tuple - x, y coordinates to convert
    :param figs: int - number of figures to output
    :return gridref: str - BNG grid reference

    Examples:

    Single value
    >>> lonlat2osgrid((-5.21469, 49.96745))

    For multiple values, use Python's zip function and list comprehension
    >>> x = [-5.21469, -5.20077, -5.18684]
    >>> y = [49.96745, 49.96783, 49.96822]
    >>> [lonlat2osgrid(coords, figs=4) for coords in zip(x, y)]
    """
    # Validate input
    bad_input_message = (
        f"Valid inputs are x, y tuple e.g. (-5.21469, 49.96783),"
        f" or list of x, y tuples. [{coords}]"
    )

    if not isinstance(coords, tuple):
        raise BNGError(bad_input_message)

    try:
        # convert to WGS84 to OSGB36 (EPSG:27700)
        # pylint: disable=E0633
        transformer = Transformer.from_crs(4326, 27700, always_xy=True)
        x_coord, y_coord = transformer.transform(coords[0], coords[1])
        # pylint: enable=E0633
    except ValueError as exc:
        raise BNGError(bad_input_message) from exc

    out_of_region_message = f"Coordinate location outside UK region: {coords}"
    if (x_coord < 0) or (y_coord < 0):
        raise BNGError(out_of_region_message)

    # Calculate region and SW corner offset

    try:
        region = _regions[int(floor(x_coord / 100000.0))][
            int(floor(y_coord / 100000.0))
        ]
        x_offset, y_offset = _offset_map[region]
    except IndexError as exc:
        raise BNGError(out_of_region_message) from exc

    # Format the output based on figs
    templates = {
        4: "{}{:02}{:02}",
        6: "{}{:03}{:03}",
        8: "{}{:04}{:04}",
        10: "{}{:05}{:05}",
    }
    factors = {4: 1000.0, 6: 100.0, 8: 10.0, 10: 1.0}
    try:  # Catch bad number of figures
        coords = templates[figs].format(
            region,
            int(floor((x_coord - x_offset) / factors[figs])),
            int(floor((y_coord - y_offset) / factors[figs])),
        )
    except KeyError as exc:
        raise BNGError("Valid inputs for figs are 4, 6, 8 or 10") from exc

    return coords


def osgrid2lonlat(gridref, epsg=None):
    """
    Convert British National Grid references to OSGB36 numeric coordinates.
    Grid references can be 4, 6, 8 or 10 figures.

    Input parameters
    ----------------
    :param gridref (str): BNG grid reference
    :param epsg (int): EPSG code

    Return:
    ------
    :coords(tuple): (x, y) coordinates

    Examples:

    Single value
    >>> osgrid2lonlat('NT2755072950', epsg=27700)
    (327550, 672950)

    For multiple values, use Python's zip function and list comprehension
    >>> gridrefs = ['HU431392', 'SJ637560', 'TV374354']
    >>> x, y = zip(*[osgrid2lonlat(g, epsg=27700) for g in gridrefs])
    >>> x
    (443100, 363700, 537400)
    >>> y
    (1139200, 356000, 35400)
    """
    # Validate input
    bad_input_message = (
        f"Valid gridref inputs are 4, 6, 8 or 10-fig references as strings "
        f'e.g. "NN123321", or lists/tuples/arrays of strings. \'[{gridref}]'
    )

    try:
        gridref = gridref.upper()
        pattern = r"^([A-Z]{2})(\d{4}|\d{6}|\d{8}|\d{10})$"
        match = re.match(pattern, gridref)
    except (TypeError, AttributeError) as exc:
        # Non-string values will throw error
        raise BNGError(bad_input_message) from exc

    if not match:
        raise BNGError(bad_input_message)

    # Extract data from gridref
    region, coords = match.groups()

    # Get offset from region
    try:
        _offset_map[region]
    except KeyError as exc:
        raise BNGError(f"Invalid 100 km grid square code: {region}") from exc

    # Get easting and northing from text and convert to coords

    easting = int(coords[: (len(coords) // 2)])
    northing = int(coords[(len(coords) // 2) :])
    scale_factor = 10 ** (5 - (len(coords) // 2))
    x_coord = int(easting * scale_factor + _offset_map[region][0])
    y_coord = int(northing * scale_factor + _offset_map[region][1])

    if epsg is None:
        return x_coord, y_coord

    try:
        transformer = Transformer.from_crs(27700, epsg, always_xy=True)
        return transformer.transform(x_coord, y_coord)
    except Exception as exc:
        raise BNGError("Invalid EPSG code provided") from exc


def nearest(item, valuelist):
    """
    Find nearest value to item in valuelist
    """
    return min(valuelist, key=lambda x: abs(x - item))


def get_dtm_values(parcel_os_code, app_config):
    """
    Query the DTM database based on longitude and latitude to retrieve
    elevation, slope and aspect data.
    The output of this function is a dictionary with the following keys:
    'x', 'y', 'elevation', 'slope', 'aspect'
    """
    # pylint: disable=R0914
    db_name = app_config.dem_parameters["db_name"]
    db_user = app_config.dem_parameters["username"]
    db_password = app_config.dem_parameters["password"]

    conn = None

    # retrieve lon, lat from parcel_os_code and create a bounding box to
    # find the closest 50m grid cell in the DEM
    lon, lat = osgrid2lonlat(parcel_os_code)
    lon_min, lon_max, lat_min, lat_max = lon - 50, lon + 50, lat - 50, lat + 50
    try:
        conn = psycopg2.connect(
            user=db_user,
            password=db_password,
            database=db_name,
            host="127.0.0.1",
            port="5432",
        )
        conn.autocommit = True
        cur = conn.cursor()
        sql = f"""
            SELECT
                terrain.x,
                terrain.y,
                terrain.val,
                terrain.slope,
                terrain.aspect
            FROM dtm.dtm_slope_aspect AS terrain
            WHERE terrain.x BETWEEN {lon_min} AND {lon_max}
            AND terrain.y BETWEEN {lat_min} AND {lat_max};
        """
        cur.execute(sql)
        sql_return = cur.fetchall()
        lon_lst = [x[0] for x in sql_return]
        lat_lst = [x[1] for x in sql_return]
        closest_lon, closest_lat = nearest(lon, lon_lst), nearest(lat, lat_lst)
        ind = [
            i for i, x in enumerate(sql_return) if x[0:2] == (closest_lon, closest_lat)
        ]
        dtm_vals = sql_return[ind[0]]
        dict_keys = ["x", "y", "elevation", "slope", "aspect"]
        dtm_dict = dtm_dict = dict(zip(dict_keys, dtm_vals))
        return dtm_dict
    except psycopg2.DatabaseError as error:
        print(error)
        return None
    finally:
        if conn is not None:
            conn.close()


def read_parcel_data(gid, folder):
    """
    Read UKCEH parcel data retrieved from a series of
    files contained in "folder".

    Parameters
    ----------
    :param gid (int): parcel ID. It must one of the
        parcel IDs of the 2021 CEH land cover map
    :param main_folder (str): path where all the json
        files containing parcel data are stored.

    Return
    ------
    rp (List): a list containing [lon, lat] corrdinates in
        EPSG:4326 of the centroid of the parcel with gid =
        inputted gid
    """
    metadata = {}
    for root, _, files in os.walk(folder):
        for filename in files:
            if filename.endswith("_meta.json"):
                filepath = os.path.join(root, filename)
                with open(filepath, encoding="utf-8") as file:
                    data = json.load(file)
                    for feature in data["features"]:
                        parcel_id = feature["gid"]
                        coords = feature["rp"]
                        metadata[parcel_id] = coords
    return metadata.get(gid, "Parcel not found")

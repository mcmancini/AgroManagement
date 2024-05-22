# -*- coding: utf-8 -*-
# Copyright (c) 2024 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), May 2024
# ===================================================
"""
LaiGenerator: a class containing a series of methods to 
generate LAI values from the comination of Sentinel 1, Sentinel 2,
and ERA 5 reanalysis data using a back-filling algorithm developed by 
Myrgiotis and Vasilis (2021) https://datashare.ed.ac.uk/handle/10283/4086.
This has been adapted and significanlty altered to make it work and have 
better formatting. It is still far from good (for example, I believe all
the procedures here implemented using subprocess.call() can be performed
using snappy); there are a number of things to check (processing of sentinel
2 does not have masking based on the GeoJSON geometry passed as input; for
some reason after backfilling the script takes the first LAI value for each
week, rather than doing for example mean resampling -- I changed that, but
check if it is the right thing to do; a ton of Pandas warnings and warnings
from SNAP are returned during execution);
There are a number of improvements to make (due to issues with my PC, I am
calling powershell in most subprocess.call(), which will not work on Mac and 
Linux).
These issues will be addressed in future updates (check issues in
https://github.com/LEEP-Modelling-Team/AgroManagement/tree/iss_24_estimate_LAI).
"""
import datetime
import glob
import json
import os
import subprocess
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterstats
import xarray as xr
from osgeo import gdal
from osgeo_utils.gdal_calc import Calc
from sentinelhub import BBoxSplitter
from shapely.geometry import shape
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

gdal.UseExceptions()
warnings.filterwarnings("ignore")
warnings.filterwarnings(
    "ignore", message="TIFFReadDirectory:Sum of Photometric type-related color channels and ExtraSamples doesn't match SamplesPerPixel.")


class LaiGenerator:
    """
    ---------------------------------------------------------------
    Production of continuous weekly LAI time-series by fusing earth
    observation data from the Sentinel-1 and Sentinel-2 systems
    ---------------------------------------------------------------

    Methods defined here:

    - s1_to_vvvh(): Processes Sentinel-1 granules into backscatter VV/VH
        intensity maps. Output shows data for the area-of-interest only
    - s2_to_lai(): Processes Sentinel-2 images/tiles into Leaf Area Index.
        Output shows data for the entire S2 tile.
    - lai_ts_creation(): Collects the outputs of S1_to_VVVH() and S2_to_LAI().
        Creates time-series of vapour pressure deficit (VPD) using data
        downloaded from ECMWF ERA-5. Returns weekly continuous LAI time-series
        for a given field and time-period.

    User inputs :

    jsonloc (str)          : path to a field's .geojson file
    filename (str)         : name of field
    workingdir (str)       : working directory (will be created if not existing)
    s1_dir (str)           : folder containing the raw Sentinel-1 data
    s2_dir (str)           : folder containing the raw Sentinel-2 data
    era_dir (str)          : folder containing ECMWF ERA-5 VPD times-series
    startdate (str)        : first day of simulations (e.g. 2018-01-01)
    enddate (str)          : last day of simulations (e.g. 2020-12-31)
    snap_graphs_dir (str)  : directory containing the ESA SNAP graphs
    snap_gtp_dir (str)     : directory containing the ESA SNAP gpt app/exe

    Notes :

    - The Sentinel-2 data have to be atmospherically corrected
        (i.e. L2A data product) and have a spatial resolution <= 20m
    - s1_dir and s2_dir should contain S2_tiles/S1_granules for the area and
        time-period of interest ONLY
    - See README.md for template .py script for obtaining VPD data from ECMWF.
        Keep the ERA5_YYYY.nc file naming convention
    - The temporal coverage of the provided S1, S2 and VPD data should correspond
        to the examined time period (i.e. startdate to enddate)
    - If S1 and S2 data are pre-processed to VV/VH and LAI the file naming convention
        should follow that used in this script (see SNAP-calling lines 247 and 367)
    --------------------------------------
    """
    # pylint: disable=R0913
    def __init__(
        self,
        jsonloc,
        filename,
        workingdir,
        s1_dir,
        s2_dir,
        era_dir,
        snap_graphs_dir,
        snap_gtp_dir,
        startdate,
        enddate,
    ):

        self.jsonloc = jsonloc
        self.filename = filename
        self.workingdir = workingdir
        self.s1_dir = s1_dir
        self.s2_dir = s2_dir
        self.era_dir = era_dir
        self.snap_graphs_dir = snap_graphs_dir
        self.snap_gtp_dir = snap_gtp_dir
        self.startdate = startdate
        self.enddate = enddate

        os.makedirs(self.workingdir, exist_ok=True)

    def s1_to_vvvh(self):
        """
        Process all available Sentinel-1 granules into VV/VH backscatter
        intensity (10m resolution)
        - Please keep the s1_dir inside the workingdir
        - Uses ESA SNAP to produce VV/VH db from Sentinel-1 data
        - For the data processing pipeline see Truckenbrodt et al 2019
            (https://doi.org/10.3390/data4030093)
        """

        fieldpolygonloc = gpd.read_file(self.jsonloc)
        poly = str(fieldpolygonloc.geometry.iloc[0])
        os.chdir(self.s1_dir)
        folds = glob.glob("*")

        for _, item in enumerate(folds):
            document = f"""\
			<graph id="Graph">
			<version>1.0</version>
			<node id="Read">
				<operator>Read</operator>
				<sources/>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<file>{self.s1_dir}/{item}</file>
					<formatName>SENTINEL-1</formatName>
				</parameters>
			</node>
			<node id="Apply-Orbit-File">
				<operator>Apply-Orbit-File</operator>
				<sources>
					<sourceProduct refid="Read"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<orbitType>Sentinel Precise (Auto Download)</orbitType>
					<polyDegree>3</polyDegree>
					<continueOnFail>false</continueOnFail>
				</parameters>
			</node>
			<node id="Remove-GRD-Border-Noise">
				<operator>Remove-GRD-Border-Noise</operator>
				<sources>
					<sourceProduct refid="Apply-Orbit-File"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<selectedPolarisations>VH,VV</selectedPolarisations>
					<borderLimit>500</borderLimit>
					<trimThreshold>0.5</trimThreshold>
				</parameters>
			</node>
			<node id="ThermalNoiseRemoval">
				<operator>ThermalNoiseRemoval</operator>
				<sources>
					<sourceProduct refid="Remove-GRD-Border-Noise"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<selectedPolarisations>VH,VV</selectedPolarisations>
					<removeThermalNoise>true</removeThermalNoise>
					<reIntroduceThermalNoise>false</reIntroduceThermalNoise>
				</parameters>
			</node>
			<node id="Calibration">
				<operator>Calibration</operator>
				<sources>
					<sourceProduct refid="ThermalNoiseRemoval"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<sourceBands/>
					<auxFile>Product Auxiliary File</auxFile>
					<externalAuxFile/>
					<outputImageInComplex>false</outputImageInComplex>
					<outputImageScaleInDb>false</outputImageScaleInDb>
					<createGammaBand>false</createGammaBand>
					<createBetaBand>false</createBetaBand>
					<selectedPolarisations>VH,VV</selectedPolarisations>
					<outputSigmaBand>true</outputSigmaBand>
					<outputGammaBand>false</outputGammaBand>
					<outputBetaBand>false</outputBetaBand>
				</parameters>
			</node>
			<node id="Multilook">
				<operator>Multilook</operator>
				<sources>
					<sourceProduct refid="Calibration"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<sourceBands>Sigma0_VH,Sigma0_VV</sourceBands>
					<nRgLooks>1</nRgLooks>
					<nAzLooks>1</nAzLooks>
					<outputIntensity>true</outputIntensity>
					<grSquarePixel>true</grSquarePixel>
				</parameters>
			</node>
			<node id="Terrain-Correction">
				<operator>Terrain-Correction</operator>
				<sources>
					<sourceProduct refid="Multilook"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<sourceBands>Sigma0_VH,Sigma0_VV</sourceBands>
					<demName>SRTM 3Sec</demName>
					<externalDEMFile/>
					<externalDEMNoDataValue>0.0</externalDEMNoDataValue>
					<externalDEMApplyEGM>true</externalDEMApplyEGM>
					<demResamplingMethod>BILINEAR_INTERPOLATION</demResamplingMethod>
					<imgResamplingMethod>BILINEAR_INTERPOLATION</imgResamplingMethod>
					<pixelSpacingInMeter>10.0</pixelSpacingInMeter>
					<pixelSpacingInDegree>8.983152841195215E-5</pixelSpacingInDegree>
					<mapProjection>GEOGCS[&quot;WGS84(DD)&quot;, 
			DATUM[&quot;WGS84&quot;, 
				SPHEROID[&quot;WGS84&quot;, 6378137.0, 298.257223563]], 
			PRIMEM[&quot;Greenwich&quot;, 0.0], 
			UNIT[&quot;degree&quot;, 0.017453292519943295], 
			AXIS[&quot;Geodetic longitude&quot;, EAST], 
			AXIS[&quot;Geodetic latitude&quot;, NORTH]]</mapProjection>
					<alignToStandardGrid>false</alignToStandardGrid>
					<standardGridOriginX>0.0</standardGridOriginX>
					<standardGridOriginY>0.0</standardGridOriginY>
					<nodataValueAtSea>true</nodataValueAtSea>
					<saveDEM>false</saveDEM>
					<saveLatLon>false</saveLatLon>
					<saveIncidenceAngleFromEllipsoid>false</saveIncidenceAngleFromEllipsoid>
					<saveLocalIncidenceAngle>false</saveLocalIncidenceAngle>
					<saveProjectedLocalIncidenceAngle>false</saveProjectedLocalIncidenceAngle>
					<saveSelectedSourceBand>true</saveSelectedSourceBand>
					<outputComplex>false</outputComplex>
					<applyRadiometricNormalization>false</applyRadiometricNormalization>
					<saveSigmaNought>false</saveSigmaNought>
					<saveGammaNought>false</saveGammaNought>
					<saveBetaNought>false</saveBetaNought>
					<incidenceAngleForSigma0>Use projected local incidence angle from DEM</incidenceAngleForSigma0>
					<incidenceAngleForGamma0>Use projected local incidence angle from DEM</incidenceAngleForGamma0>
					<auxFile>Latest Auxiliary File</auxFile>
					<externalAuxFile/>
				</parameters>
			</node>
			<node id="Speckle-Filter">
				<operator>Speckle-Filter</operator>
				<sources>
					<sourceProduct refid="Terrain-Correction"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<sourceBands>Sigma0_VH,Sigma0_VV</sourceBands>
					<filter>Lee Sigma</filter>
					<filterSizeX>3</filterSizeX>
					<filterSizeY>3</filterSizeY>
					<dampingFactor>2</dampingFactor>
					<estimateENL>true</estimateENL>
					<enl>1.0</enl>
					<numLooksStr>1</numLooksStr>
					<windowSize>7x7</windowSize>
					<targetWindowSizeStr>3x3</targetWindowSizeStr>
					<sigmaStr>0.9</sigmaStr>
					<anSize>50</anSize>
				</parameters>
			</node>
			<node id="Subset">
				<operator>Subset</operator>
				<sources>
					<sourceProduct refid="Speckle-Filter"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<sourceBands>Sigma0_VH,Sigma0_VV</sourceBands>
					<region>0,0,0,0</region>
					<referenceBand/>
					<geoRegion>{poly}</geoRegion>
					<subSamplingX>1</subSamplingX>
					<subSamplingY>1</subSamplingY>
					<fullSwath>false</fullSwath>
					<tiePointGridNames/>
					<copyMetadata>true</copyMetadata>
				</parameters>
			</node>
			<node id="Write">
				<operator>Write</operator>
				<sources>
					<sourceProduct refid="Subset"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<file>{self.s1_dir}/processed/Subset_{item[:-4]}_Orb_NR_Cal_ML_TC_dB.tif</file>
					<formatName>GeoTIFF</formatName>
				</parameters>
			</node>
			<applicationData id="Presentation">
				<Description/>
				<node id="Read">
					<displayPosition x="37.0" y="134.0"/>
				</node>
				<node id="Apply-Orbit-File">
					<displayPosition x="157.0" y="75.0"/>
				</node>
				<node id="Remove-GRD-Border-Noise">
					<displayPosition x="275.0" y="161.0"/>
				</node>
				<node id="ThermalNoiseRemoval">
					<displayPosition x="475.0" y="90.0"/>
				</node>
				<node id="Calibration">
					<displayPosition x="634.0" y="166.0"/>
				</node>
				<node id="Multilook">
					<displayPosition x="734.0" y="92.0"/>
				</node>
				<node id="Terrain-Correction">
					<displayPosition x="854.0" y="168.0"/>
				</node>
				<node id="Speckle-Filter">
					<displayPosition x="981.0" y="112.0"/>
				</node>
				<node id="Subset">
					<displayPosition x="1101.0" y="136.0"/>
				</node>
				<node id="Write">
					<displayPosition x="1195.0" y="181.0"/>
				</node>
			</applicationData>
			</graph>
			"""

            os.chdir(self.snap_graphs_dir)
            with open("S1_to_VVVH.xml", "w", encoding="utf-8") as file1:
                file1.write(document)

            xml_file_path = f"{self.snap_graphs_dir}S1_to_VVVH.xml"
            subprocess.call([f"{self.snap_gtp_dir}gpt.exe", xml_file_path])

    def s2_to_lai(self):
        """
        --------------------------------------
        Process all available Sentinel-2 images/tiles into maps of Leaf
        Area Index (LAI)
        --------------------------------------
        - Please keep the s2_dir inside the workingdir
        - This functions performs the following processes :
        --- Apply ESA SNAP resampling and biophysical calculator to produce LAI data
        --- Reproject the resulting LAI maps to EPSG:4326
        --- Remove cloud-affected pixels from final .tif
        *** Final .tif has _p2 appended to its name
        *** The final LAI maps' spatial resolution can be adjusted by changing the
            <targetResolution> (line 332) -- currently 40m
        --------------------------------------
        """

        # fieldpolygonloc = gpd.read_file(self.jsonloc)
        # poly = str(fieldpolygonloc.geometry.iloc[0])
        os.chdir(self.s2_dir)
        subprocess.call("rm -r processed", shell=True)
        folds = glob.glob("*")

        for _, item in enumerate(folds):
            document = f"""\
            <graph id="Graph">
            <version>1.0</version>
            <node id="Read">
                <operator>Read</operator>
                <sources/>
                <parameters class="com.bc.ceres.binding.dom.XppDomElement">
                <useAdvancedOptions>false</useAdvancedOptions>
                <file>{self.s2_dir}{item}</file>
                <copyMetadata>true</copyMetadata>
                <bandNames/>
                <pixelRegion>0,0,10980,10980</pixelRegion>
                <maskNames/>
                </parameters>
            </node>
            <node id="Resample">
                <operator>Resample</operator>
                <sources>
                <sourceProduct refid="Read"/>
                </sources>
                <parameters class="com.bc.ceres.binding.dom.XppDomElement">
                <referenceBand/>
                <targetWidth/>
                <targetHeight/>
                <targetResolution>40</targetResolution>
                <upsampling>Nearest</upsampling>
                <downsampling>Mean</downsampling>
                <flagDownsampling>First</flagDownsampling>
                <resamplingPreset/>
                <bandResamplings/>
                <resampleOnPyramidLevels>true</resampleOnPyramidLevels>
                </parameters>
            </node>
            <node id="BiophysicalOp">
                <operator>BiophysicalOp</operator>
                <sources>
                <sourceProduct refid="Resample"/>
                </sources>
                <parameters class="com.bc.ceres.binding.dom.XppDomElement">
                <sensor>S2A</sensor>
                <computeLAI>true</computeLAI>
                <computeFapar>false</computeFapar>
                <computeFcover>false</computeFcover>
                <computeCab>false</computeCab>
                <computeCw>false</computeCw>
                </parameters>
            </node>
            <node id="Write">
				<operator>Write</operator>
				<sources>
				  <sourceProduct refid="BiophysicalOp"/>
				</sources>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
				  <file>{self.s2_dir}/processed/{item[:-4]}.tif</file>
				  <formatName>GeoTIFF</formatName>
				</parameters>
            </node>
            <applicationData id="Presentation">
				<Description/>
				<node id="Read">
						<displayPosition x="37.0" y="134.0"/>
				</node>
				<node id="BiophysicalOp">
				  <displayPosition x="308.0" y="102.0"/>
				</node>
				<node id="Resample">
				  <displayPosition x="184.0" y="98.0"/>
				</node>
				<node id="Write">
						<displayPosition x="455.0" y="135.0"/>
				</node>
            </applicationData>
            </graph>
            """

            os.chdir(self.snap_graphs_dir)
            with open("S2_to_LAI.xml", "w", encoding="utf-8") as file1:
                file1.write(document)
            xml_file_path = f"{self.snap_graphs_dir}S2_to_LAI.xml"
            subprocess.call([f"{self.snap_gtp_dir}gpt.exe", xml_file_path])
            try:
                directory = (
                    f"{self.snap_graphs_dir[:-7]}var/" f"cache/s2tbx/l2a-reader/10.0.1"
                )
                powershell_command = (
                    f"Get-ChildItem -Path {directory} " f"| Remove-Item -Recurse -Force"
                )

                subprocess.call(["powershell", powershell_command])
            except subprocess.CalledProcessError as e:
                print("Error:", e)

        # Reproject LAI tiffs and remove cloud pixels
        os.chdir(f"{self.s2_dir}processed")
        folders = []
        for file in glob.glob("*.tif"):
            folders.append(file)
        folders.sort()
        for _, item in enumerate(folders):
            subprocess.call(
                f"gdalwarp {item} {item[:-4]}_p1.tif "
                f"-s_srs EPSG:32630 -t_srs EPSG:4326",
            )

            subprocess.call(["powershell", "-Command", f"Remove-Item -Path '{item}'"])

            Calc(
                calc="A*(B==0)",
                A=f"{item[:-4]}_p1.tif",
                A_band=1,
                B=f"{item[:-4]}_p1.tif",
                B_band=2,
                outfile=f"{item[:-4]}_p2.tif",
                NoDataValue=0,
            )

            subprocess.call(["powershell", "-Command", "Remove-Item -Path ./*s_p1.tif"])

    # pylint: disable=R0912, R0915
    def lai_ts_creation(self):
        """
        ---------------------------------------------------------------------------
        Collect the processed S1-VV/VH and S2-LAI data,
        calculate VPD and produce continuous weekly field-scale LAI (m2/m2)
        ---------------------------------------------------------------------------
        uses: (1) S2-based LAI, (2) S1-based VV and VH backscatter
            intensity data (3) ECMWF ERA5-based temperature and dewpoint data
        returns: a numpy array of weekly field-mean LAI (m2.m-2) for the
            given period and a numpy array with the R2 of the Random Forest model
        ---------------------------------------------------------------------------
        """

        fieldpolygonloc = self.jsonloc

        ### S2-LAI data directory
        os.chdir(f"{self.s2_dir}processed")
        folders_s2 = []
        for file in glob.glob("*_p2.tif"):
            folders_s2.append(file)
        folders_s2.sort()

        ### Split fields in sub-fields
        with open(fieldpolygonloc, encoding="utf-8") as f:
            js = json.load(f)
        for feature in js["features"]:
            polygon = shape(feature["geometry"])
        bbox_splitter = BBoxSplitter(
            [polygon], 4326, (5, 5)
        )  # bounding box will be split into x-times-x bounding boxes

        ### collect S2 LAI data per sub-field
        sentinel_2_df = pd.DataFrame(columns=["box", "date", "lai", "lai_std"])
        for y, _ in enumerate(folders_s2):
            for i, _ in enumerate(bbox_splitter.get_bbox_list()):
                listofzones_lai = rasterstats.zonal_stats(
                    (bbox_splitter.bbox_list[i]).geometry.wkt,
                    folders_s2[y],
                    stats=["mean", "std", "count"],
                    band=1,
                    nodata=0,
                )
                new_row = pd.DataFrame(
                    [
                        {
                            "box": int(i),
                            "date": datetime.datetime.strptime(
                                folders_s2[y][11:19], "%Y%m%d"
                            ),
                            "lai": listofzones_lai[0]["mean"],
                            "lai_std": listofzones_lai[0]["std"],
                        }
                    ]
                )
                sentinel_2_df = pd.concat([sentinel_2_df, new_row], ignore_index=True)

        sentinel_2_df.index = sentinel_2_df.date
        sentinel_2_df = sentinel_2_df.sort_index()

        ### S1 SAR data directory
        os.chdir(f"{self.s1_dir}processed")
        folders_s1 = []
        for file in glob.glob("*.tif"):
            folders_s1.append(file)
        folders_s1.sort()

        ### Collect S1 backscatter data per sub-field
        sentinel_1_df = pd.DataFrame(
            columns=["box", "date", "band1", "band1_std", "band2", "band2_std"]
        )
        for y, folder in enumerate(folders_s1):
            for ii, bbox in enumerate(bbox_splitter.get_bbox_list()):
                band1 = rasterstats.zonal_stats(
                    bbox.geometry.wkt,
                    f"{self.s1_dir}processed/{folder}",
                    stats=["mean", "std"],
                    band=1,
                    nodata=0,
                )

                band2 = rasterstats.zonal_stats(
                    bbox.geometry.wkt,
                    folder,
                    stats=["mean", "std"],
                    band=2,
                    nodata=0,
                )

                new_row = pd.DataFrame(
                    [
                        {
                            "box": int(ii),
                            "date": datetime.datetime.strptime(folder[24:32], "%Y%m%d"),
                            "band1": band1[0]["mean"],
                            "band1_std": band1[0]["std"],
                            "band2": band2[0]["mean"],
                            "band2_std": band2[0]["std"],
                        }
                    ]
                )
                sentinel_1_df = pd.concat([sentinel_1_df, new_row], ignore_index=True)

        sentinel_1_df.index = sentinel_1_df.date
        sentinel_1_df = sentinel_1_df.sort_index()
        sentinel_1_df["bandratio"] = sentinel_1_df.band1 / sentinel_1_df.band2

        ### Merge S1 and S2 per sub-field
        merged_sentinel = sentinel_1_df["2017":"2019"]
        merged_sentinel = merged_sentinel.sort_index()
        merged_sentinel["lai"] = np.nan
        merged_sentinel["lai_std"] = np.nan
        sentinel_2_v2 = sentinel_2_df[sentinel_2_df.lai > 0]
        merged_sentinel = merged_sentinel[
            merged_sentinel.date.isin(sentinel_2_v2.date.unique())
        ]
        for ii in range(len(set(merged_sentinel.date))):
            for b in range(len(set(sentinel_2_v2.box))):
                date_mask = merged_sentinel.date == merged_sentinel.date.unique()[ii]
                box_mask = merged_sentinel.box == sentinel_2_v2.box.unique()[b]

                df_filter = (
                    sentinel_2_v2.date == merged_sentinel.date.unique()[ii]
                ) & (sentinel_2_v2.box == sentinel_2_v2.box.unique()[b])
                filtered_s2_v2 = sentinel_2_v2[df_filter]

                if not filtered_s2_v2.empty:
                    merged_sentinel.loc[date_mask & box_mask, "lai"] = filtered_s2_v2[
                        "lai"
                    ].values[0]

                    merged_sentinel.loc[date_mask & box_mask, "lai_std"] = (
                        filtered_s2_v2["lai_std"].values[0]
                    )
        del b, ii

        ### Load and process met data
        os.chdir(f"{self.workingdir}/{self.era_dir}")
        year_first = self.startdate[:4]
        year_last = self.enddate[:4]
        list_of_years = np.arange(int(year_first), int(year_last) + 1)
        folders_met = []
        for file in glob.glob("*.nc"):
            if file[5:9] in str(list_of_years):
                folders_met.append(file)
        folders_met.sort()

        ### load field location info
        shapefile = gpd.read_file(fieldpolygonloc)
        shapefile["lat"] = shapefile.geometry.centroid.y.iloc[0]
        shapefile["lon"] = shapefile.geometry.centroid.x.iloc[0]

        # fill a era_df pandas dataframe with the climate data needed to calculate
        # VPD - resample from hourly to daily
        era_df = pd.DataFrame(columns=["date", "vpd"])
        for ii, folder_met in enumerate(folders_met):
            ds = xr.open_dataset(folder_met)
            dsloc = ds.sel(
                longitude=float(shapefile["lon"]),  # polygon centroid lon
                latitude=float(shapefile["lat"]),  # polygon centroid lat
                method="nearest",
            )
            df = dsloc.to_dataframe()

            # Unit conversion
            df["t2m"] -= 273.15
            df["d2m"] -= 273.15

            # Relative humidity calculation
            df["RH"] = 100 * (
                np.exp((17.625 * df["d2m"]) / (243.04 + df["d2m"]))
                / np.exp((17.625 * df["t2m"]) / (243.04 + df["t2m"]))
            )

            # Vapor pressure deficit calculation
            df["VPD"] = (1 - (df["RH"] / 100)) * (
                610.7 * 10 ** (7.5 * df["t2m"] / (237.3 + df["t2m"]))
            )

            df = df.reset_index()
            df.index = pd.date_range(
                f"{folders_met[ii][5:9]}-01-01", periods=len(df), freq="h"
            )

            # Resample from hourly to daily
            met_data = []

            for t in range(int(len(df) / float(24))):
                met_data.append(
                    {
                        "date": pd.date_range(
                            f"{folders_met[ii][5:9]}-01-01",
                            periods=int(len(df) / float(8)),
                            freq="D",
                        )[t],
                        "vpd": df["VPD"].resample("D", label="left").mean().iloc[t],
                    }
                )
            del t
            era_df = pd.concat([era_df, pd.DataFrame(met_data)], ignore_index=True)

            era_df.index = era_df["date"]
            era_df["DOY"] = era_df.index.dayofyear

        ### Add met info to S1+S2 dataframe
        merged_sentinel = merged_sentinel.dropna()
        merged_sentinel["DOY"] = merged_sentinel.index.dayofyear
        merged_sentinel["vpd"] = era_df.vpd

        ### Train Random Forest algorithm using box/subfield data
        x_train, x_test, y_train, y_test = train_test_split(
            merged_sentinel[["band1", "band2", "DOY", "vpd"]],
            merged_sentinel.lai,
            test_size=0.2,
            random_state=0,
        )
        rf = RandomForestRegressor(n_estimators=100)
        rf.fit(x_train, y_train)
        rf_score = rf.score(x_test, y_test)
        np.save(
            f"{self.workingdir}{self.filename}_RF_R2={rf_score}.npy",
            np.array([rf_score]),
        )

        ### Fill S1 dataframe with RF predcited LAI
        sentinel_1_df = sentinel_1_df[self.startdate :]
        sentinel_1_df["DOY"] = sentinel_1_df.index.dayofyear
        # daily average cross all boxes
        sentinel_1_df = sentinel_1_df.resample("D").median()
        sentinel_1_df = sentinel_1_df.dropna()
        sentinel_1_df["vpd"] = era_df.vpd
        sentinel_1_df["rf_LAI"] = np.nan
        sentinel_1_df["rf_LAI"] = rf.predict(
            sentinel_1_df[["band1", "band2", "DOY", "vpd"]]
        )

        ## Use mean field VV/VH to RF-predict mean field LAI
        daily_rf_lai = pd.DataFrame(
            columns=["rf_LAI"],
            index=pd.date_range(start=self.startdate, end=self.enddate, freq="D"),
        )
        daily_rf_lai["rf_LAI"] = round(sentinel_1_df["rf_LAI"], 2)
        daily_rf_lai.loc[0, "rf_LAI"] = 0
        daily_rf_lai = daily_rf_lai.interpolate(
            "linear"
        )  # interpolated RF LAI time series
        daily_rf_lai.index = pd.to_datetime(daily_rf_lai.index)
        daily_rf_lai = daily_rf_lai.iloc[:-1]
        lai_timeseries = daily_rf_lai.rf_LAI.resample("7D", label="right").mean()[:-1]
        lai_timeseries_df = lai_timeseries.to_frame().reset_index()
        lai_timeseries_df.rename(columns={"index": "date"})
        lai_timeseries_df.to_csv(
            f"{self.workingdir}/lai_timeseries_{self.filename}.csv", index=False
        )


# pylint: enable=R0912, R0913, R0915
warnings.filterwarnings("default")

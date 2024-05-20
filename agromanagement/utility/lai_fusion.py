import json
import pandas as pd
import numpy as np
from shapely.geometry import shape, Point
import datetime
import glob
import matplotlib.pyplot as plt
import os
import subprocess
import rasterstats
import geopandas as gpd
import pandas as pd
from osgeo import gdal
from osgeo_utils.gdal_calc import Calc
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sentinelhub import BBoxSplitter

gdal.UseExceptions()

class LAI_S1S2_fusion:
    """
    --------------------------------------
    Production of continuous weekly LAI time-series by fusing earth observation data from the Sentinel-1 and Sentinel-2 systems
    --------------------------------------

    Description of functions :

    - S1_to_VVVH() : Processes Sentinel-1 granules into backscatter VV/VH intensity maps. Output shows data for the area-of-interest only
    - S2_to_LAI()  : Processes Sentinel-2 images/tiles into Leaf Area Index. Output shows data for the entire S2 tile.
    - lai_ts_creation() : Collects the outputs of S1_to_VVVH() and S2_to_LAI(). Creates time-series of vapour pressure deficit (VPD)
    using data downloaded from ECMWF ERA-5. Returns weekly continuous LAI time-series for a given field and time-period.

    User inputs :

    jsonloc (str)         : path to a field's .geojson file
    FNAME (str)           : name of field
    workingdir (str)      : working directory (will be created if not existing)
    s1dir (str)           : folder containing the raw Sentinel-1 data (should be located inside workingdir)
    s2dir (str)           : folder containing the raw Sentinel-2 data (should be located inside workingdir)
    ECMWF_ERA5_VPD (str)  : folder containing ECMWF ERA-5 VPD times-series (should be located inside workingdir)
    startdate (str)       : first day of simulations (e.g. 2018-01-01)
    enddate (str)         : last day of simulations (e.g. 2020-12-31)
    snap_graphs_dir (str) : directory containing the ESA SNAP graphs
    snap_gtp_dir (str)    : directory containing the ESA SNAP gpt app/exe

    Notes :

    - The Sentinel-2 data have to be atmospherically corrected (i.e. L2A data product) and have a spatial resolution <= 20m
    - s1dir and s2dir should contain S2_tiles/S1_granules for the area and time-period of interest ONLY
    - See README.md for template .py script for obtaining VPD data from ECMWF. Keep the ERA5_YYYY.nc file naming convention
    - The temporal coverage of the provided S1, S2 and VPD data should correspond to the examined time period (i.e. startdate to enddate)
    - If S1 and S2 data are pre-processed to VV/VH and LAI the file naming convention should follow that used in this script (see SNAP-calling lines 247 and 367)
    --------------------------------------
    """

    def __init__(
        self,
        jsonloc,
        FNAME,
        workingdir,
        s1dir,
        s2dir,
        ECMWF_ERA5_VPD,
        snap_graphs_dir,
        snap_gtp_dir,
        startdate,
        enddate,
    ):

        self.jsonloc = jsonloc
        self.FNAME = FNAME
        self.workingdir = workingdir
        self.s1dir = s1dir
        self.s2dir = s2dir
        self.ECMWF_ERA5_VPD = ECMWF_ERA5_VPD
        self.snap_graphs_dir = snap_graphs_dir
        self.snap_gtp_dir = snap_gtp_dir
        self.startdate = startdate
        self.enddate = enddate

        os.makedirs(self.workingdir, exist_ok=True)

    def S1_to_VVVH(self):
        """
        --------------------------------------
                Process all available Sentinel-1 granules into VV/VH backscatter intensity (10m resolution)
        --------------------------------------
                - Please keep the s1dir inside the workingdir
                - Uses ESA SNAP to produce VV/VH db from Sentinel-1 data
                - For the data processing pipeline see Truckenbrodt et al 2019 (https://doi.org/10.3390/data4030093)
        --------------------------------------
        """

        fieldpolygonloc = gpd.read_file(self.jsonloc)
        poly = str(fieldpolygonloc.geometry.iloc[0])
        os.chdir(self.s1dir)
        folds = glob.glob("*")

        for _, item in enumerate(folds):
            document = f"""\
			<graph id="Graph">
			<version>1.0</version>
			<node id="Read">
				<operator>Read</operator>
				<sources/>
				<parameters class="com.bc.ceres.binding.dom.XppDomElement">
					<file>{self.s1dir}/{item}</file>
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
					<file>{self.s1dir}/processed/Subset_{item[:-4]}_Orb_NR_Cal_ML_TC_dB.tif</file>
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
            file1 = open("S1_to_VVVH.xml", "w", encoding="utf-8")
            file1.write(document)
            file1.close()

            xml_file_path = f"{self.snap_graphs_dir}S1_to_VVVH.xml"
            subprocess.call([f"{self.snap_gtp_dir}gpt.exe", xml_file_path])

    def S2_to_LAI(self):
        """
        --------------------------------------
        Process all available Sentinel-2 images/tiles into maps of Leaf Area Index (LAI)
        --------------------------------------
        - Please keep the s2dir inside the workingdir
        - This functions performs the following processes :
        --- Apply ESA SNAP resampling and biophysical calculator to produce LAI data
        --- Reproject the resulting LAI maps to EPSG:4326
        --- Remove cloud-affected pixels from final .tif
        *** Final .tif has _p2 appended to its name
        *** The final LAI maps' spatial resolution can be adjusted by changing the <targetResolution> (line 332) -- currently 40m
        --------------------------------------
        """

        fieldpolygonloc = gpd.read_file(self.jsonloc)
        poly = str(fieldpolygonloc.geometry.iloc[0])
        os.chdir(self.s2dir)
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
                <file>{self.s2dir}{item}</file>
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
				  <file>{self.s2dir}/processed/{item[:-4]}.tif</file>
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
                    f"{self.snap_graphs_dir[:-7]}var/"
                    f"cache/s2tbx/l2a-reader/10.0.1"
                )
                powershell_command = (
                    f"Get-ChildItem -Path {directory} "
                    f"| Remove-Item -Recurse -Force"
                )

                subprocess.call(["powershell", powershell_command])
            except subprocess.CalledProcessError as e:
                print("Error:", e)

        # Reproject LAI tiffs and remove cloud pixels
        os.chdir(f"{self.s2dir}processed")
        folders = []
        for file in glob.glob("*.tif"):
            folders.append(file)
        folders.sort()
        for _, item in enumerate(folders):
            subprocess.call(
                f"gdalwarp {item} {item[:-4]}_p1.tif "
                f"-s_srs EPSG:32630 -t_srs EPSG:4326",
            )

            subprocess.call([
                "powershell",
                "-Command",
                f"Remove-Item -Path '{item}'"
            ])

            Calc(
                calc="A*(B==0)",
                A=f"{item[:-4]}_p1.tif",
                A_band=1,
                B=f"{item[:-4]}_p1.tif",
                B_band=2,
                outfile=f"{item[:-4]}_p2.tif",
                NoDataValue=0
            )

            subprocess.call([
                "powershell",
                "-Command",
                "Remove-Item -Path ./*s_p1.tif"
            ])

    def lai_ts_creation(self):
        """
        --------------------------------------
        Collect the processed S1-VV/VH and S2-LAI data, calculate VPD and produce continuous weekly field-scale LAI (m2/m2)
        --------------------------------------
        uses : (1) S2-based LAI, (2) S1-based VV and VH backscatter intensity data (3) ECMWF ERA5-based temperature and dewpoint data
        returns : a numpy array of weekly field-mean LAI (m2.m-2) for the given period and a numpy array with the R2 of the Random Forest model
        --------------------------------------
        """

        fieldpolygonloc = self.jsonloc

        ### S2-LAI data directory
        os.chdir("%s/processed" % (self.s2dir))
        folders_S2 = []
        for file in glob.glob("*_p2.tif"):
            folders_S2.append(file)
        folders_S2.sort()

        ### Split fields in sub-fields
        with open(fieldpolygonloc) as f:
            js = json.load(f)
        for feature in js["features"]:
            polygon = shape(feature["geometry"])
        bbox_splitter = BBoxSplitter(
            [polygon], 4326, (5, 5)
        )  # bounding box will be split into x-times-x bounding boxes

        ### collect S2 LAI data per sub-field
        S2_DF = pd.DataFrame()
        for y in range(len(folders_S2)):
            for i in range(len((bbox_splitter.get_bbox_list()))):
                listofzones_lai = rasterstats.zonal_stats(
                    (bbox_splitter.bbox_list[i]).geometry.to_wkt(),
                    folders_S2[y],
                    stats=["mean", "std", "count"],
                    band=1,
                    nodata=0,
                )
                S2_DF = S2_DF.append(
                    {
                        "box": int(i),
                        "date": datetime.datetime.strptime(
                            folders_S2[y][19:27], "%Y%m%d"
                        ),
                        "lai": listofzones_lai[0]["mean"],
                        "lai_std": listofzones_lai[0]["std"],
                    },
                    ignore_index=True,
                )
        del y

        S2_DF.index = S2_DF.date
        S2_DF = S2_DF.sort_index()

        ### S1 SAR data directory
        os.chdir("%s/processed" % (self.s1dir))
        folders_S1 = []
        for file in glob.glob("*.tif"):
            folders_S1.append(file)
        folders_S1.sort()

        ### Collect S1 backscatter data per sub-field
        S1_DF = pd.DataFrame()
        for y in range(len(folders_S1)):
            for ii in range(len((bbox_splitter.get_bbox_list()))):
                band1 = rasterstats.zonal_stats(
                    (bbox_splitter.bbox_list[ii]).geometry.to_wkt(),
                    folders_S1[y],
                    stats=["mean", "std"],
                    band=1,
                    nodata=0,
                )
                band2 = rasterstats.zonal_stats(
                    (bbox_splitter.bbox_list[ii]).geometry.to_wkt(),
                    folders_S1[y],
                    stats=["mean", "std"],
                    band=2,
                    nodata=0,
                )
                S1_DF = S1_DF.append(
                    {
                        "box": int(ii),
                        "date": datetime.datetime.strptime(
                            folders_S1[y][24:32], "%Y%m%d"
                        ),
                        "band1": band1[0]["mean"],
                        "band1_std": band1[0]["std"],
                        "band2": band2[0]["mean"],
                        "band2_std": band2[0]["std"],
                    },
                    ignore_index=True,
                )
        del y, ii

        S1_DF.index = S1_DF.date
        S1_DF = S1_DF.sort_index()
        S1_DF["bandratio"] = S1_DF.band1 / S1_DF.band2

        ### Merge S1 and S2 per sub-field
        DF = S1_DF["2017":"2019"]
        DF = DF.sort_index()
        DF["lai"] = np.nan
        DF["lai_std"] = np.nan
        S2_DF_v2 = S2_DF[S2_DF.lai > 0]
        DF = DF[DF.date.isin(list(set(S2_DF_v2.date)))]
        for ii in range(len(set(DF.date))):
            for b in range(len(set(S2_DF_v2.box))):
                if (
                    len(
                        S2_DF_v2[
                            (S2_DF_v2.date == list(set(DF.date))[ii])
                            & (S2_DF_v2.box == list(set(S2_DF_v2.box))[b])
                        ]
                    )
                    > 0
                ):
                    DF["lai"][
                        (DF.date == list(set(DF.date))[ii])
                        & (DF.box == list(set(S2_DF_v2.box))[b])
                    ] = float(
                        S2_DF_v2.lai[
                            (S2_DF_v2.date == list(set(DF.date))[ii])
                            & (S2_DF_v2.box == list(set(S2_DF_v2.box))[b])
                        ]
                    )
                    DF["lai_std"][
                        (DF.date == list(set(DF.date))[ii])
                        & (DF.box == list(set(S2_DF_v2.box))[b])
                    ] = float(
                        S2_DF_v2.lai_std[
                            (S2_DF_v2.date == list(set(DF.date))[ii])
                            & (S2_DF_v2.box == list(set(S2_DF_v2.box))[b])
                        ]
                    )
        del b, ii

        ### Load and process met data
        os.chdir("%s/" % (self.ECMWF_ERA5_VPD))
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

        ### fill a met_DF pandas dataframe with the climate data needed to calculate VPD - resample from hourly to daily
        met_DF = pd.DataFrame(columns=["date", "vpd"])
        for ii in range(len(folders_met)):
            ds = xr.open_dataset(folders_met[ii])
            dsloc = ds.sel(
                longitude=float(shapefile["lon"]),  # polygon centroid lon
                latitude=float(shapefile["lat"]),  # polygon centroid lat
                method="nearest",
            )
            df = dsloc.to_dataframe()
            ## unit conversion
            df.t2m = df.t2m - 273.15
            df.d2m = df.d2m - 273.15
            ## relative humidity (http:/andrew.rsmas.miami.edu/bmcnoldy/Humidity.html)
            df["RH"] = 100 * (
                np.exp((17.625 * df.d2m) / (243.04 + df.d2m))
                / np.exp((17.625 * df.t2m) / (243.04 + df.t2m))
            )
            ## vapor pressure deficit (http:/cronklab.wikidot.com/calculation-of-vapour-pressure-deficit)
            df["VPD"] = (1 - (df.RH / 100)) * (
                610.7 * 10 ** (7.5 * df.t2m / (237.3 + df.t2m))
            )

            df = df.reset_index()
            df.index = pd.date_range(
                "%s-01-01" % folders_met[ii][5:9], periods=len(df), freq="H"
            )

            ### resample from hourly to daily
            for t in range(int(len(df) / float(24))):
                met_DF = met_DF.append(
                    {
                        "date": pd.date_range(
                            "%s-01-01" % folders_met[ii][5:9],
                            periods=int(len(df) / float(8)),
                            freq="D",
                        )[t],
                        "vpd": (df.VPD.resample("D", label="left").mean()).iloc[t],
                    },
                    ignore_index=True,
                )
            del t

            met_DF.index = met_DF.date
            met_DF["DOY"] = met_DF.index.dayofyear
        del ii

        ### Add met info to S1+S2 dataframe
        DF = DF.dropna()
        DF["DOY"] = DF.index.dayofyear
        DF["vpd"] = met_DF.vpd

        ### Train Random Forest algorithm using box/subfield data
        X_train, X_test, y_train, y_test = train_test_split(
            DF[["band1", "band2", "DOY", "vpd"]], DF.lai, test_size=0.2, random_state=0
        )
        rf = RandomForestRegressor(n_estimators=100)
        rf.fit(X_train, y_train)
        RF_score = rf.score(X_test, y_test)
        np.save(
            "%s/%s_RF_R2=%s.npy" % (self.workingdir, self.FNAME, RF_score),
            np.array([RF_score]),
        )  # empty .npy with the RF model's R2 as its name

        ### Fill S1 dataframe with RF predcited LAI
        S1_DF = S1_DF[self.stardate :]
        S1_DF["DOY"] = S1_DF.index.dayofyear
        S1_DF = S1_DF.resample("D").median()  # daily average cross all boxes
        S1_DF = S1_DF.dropna()
        S1_DF["vpd"] = met_DF.vpd
        S1_DF["rf_LAI"] = np.nan
        S1_DF["rf_LAI"] = rf.predict(S1_DF[["band1", "band2", "DOY", "vpd"]])

        ## Use mean field VV/VH to RF-predict mean field LAI
        daily_rfLAI = pd.DataFrame(
            columns={"emptycol"},
            index=pd.date_range(start=self.startdate, end=self.enddate, freq="D"),
        )
        daily_rfLAI["rf_LAI"] = round(S1_DF["rf_LAI"], 2)
        daily_rfLAI["rf_LAI"].iloc[0] = 0
        daily_rfLAI = daily_rfLAI.interpolate(
            "linear"
        )  # interpolated RF LAI time series
        LAI_timeseries = daily_rfLAI.rf_LAI.resample("7D", label="right").first()[
            :-1
        ]  # the LAI value returned is the first for each week within the examined period

        # Save the RF predicted weekly LAI time-series
        np.save(
            ("%s/LAI_timeseries_%s.npy" % (self.workingdir, self.FNAME)), LAI_timeseries
        )

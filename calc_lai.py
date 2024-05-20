from agromanagement.utility.lai_fusion import LAI_S1S2_fusion
from agromanagement.utility.paths import ROOT_DIR

PARCEL_PATH = "D:/Documents/GitHub/AgroManagement/resources/lcm2021_tile_11_1014_647.geojson"

SENTINEL_1_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_1/Raw/"
SENTINEL_2_UNZIPPED_FOLDER = "D:/Documents/Data/Sentinel/Sentinel_2/Raw/"
ERA_DIR = "resources/era_5/"

START_DATE = "2019-04-01"
END_DATE = "2019-04-30"

lai = LAI_S1S2_fusion(
    jsonloc=PARCEL_PATH,
    FNAME="prova",
    workingdir=ROOT_DIR,
    s1dir=SENTINEL_1_UNZIPPED_FOLDER,
    s2dir=SENTINEL_2_UNZIPPED_FOLDER,
    ECMWF_ERA5_VPD=ERA_DIR,
    snap_graphs_dir="C:/Users/mcm216/.snap/graphs/",
    snap_gtp_dir="C:/Program Files/esa-snap/bin/",
    startdate=START_DATE,
    enddate=END_DATE
)

# lai.S1_to_VVVH()
lai.S2_to_LAI()
lai.lai_ts_creation()

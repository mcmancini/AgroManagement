# A model of farm management in the UK

## Authors, afifliation and licence
**Authors**: 
- Mattia C. Mancini (m.c.mancini@exeter.ac.uk)
- Brett Day
- Amy Binner
- Muhammad Hasan
- Chris Lee
- Danny Williamson. 

**Affiliation**: Land, Environment, Economics and Policy Institute, University of Exeter.

**License**: this work is licences under...

## Content
- [Model overview](#model-overview)
    - [Model parameters](#model-parameters)
    - [Modelling crop yields](#modelling-crop-yields)
    - [Management model](#management-model)
    - [Model calibration - fitting](#yield-managment-model-calibrationfitting)
- [Setup](#setup)
    - [Pre-requisites](#pre-requisites)
    - [Installation](#installation)
- [Usage](#usage)
- [Details](#details)
- [Modules](#modules)
- [Types](#types)

## Model overview
With this project we aim to build a structural, biophysically constrained, economic model of farmers' behaviour to predict farm management practices with focus on crop planting and fertilisation in the UK, fitted and calibrated using a variety of data primarily from [DEFRA](https://www.gov.uk/government/organisations/department-for-environment-food-rural-affairs) and Earth obeservations, using a [UK implementation](https://github.com/mcmancini/UkWofost) of the [WOFOST Crop yield model](https://www.wur.nl/en/research-results/research-institutes/environmental-research/facilities-tools/software-models-and-databases/wofost.htm) from the University of Wageningen in order to estimate crop yields.

### Model parameters
The main parameters of the model are the following:
- Soil functions describing water conductivity and water retention in the soil - $\boldsymbol{\alpha}$
- Plant physiology - $\boldsymbol{\beta}$
- Weather, soil type - $\boldsymbol{w}$
- Initial soil conditions (nutrients and water content) - $\boldsymbol{s}$
- Managment choices - $\boldsymbol{x, q}$

Of these, $\boldsymbol{\alpha}$ and $\boldsymbol{\beta}$ are known but uncertain for the UK, $\boldsymbol{w}$ and $\boldsymbol{s}$ are known, and the management choices $\boldsymbol{x}$ and $\boldsymbol{q}$ are endpogenous.
### Modelling crop yields
The management model relies on the UK implementation of the WOFOST crop yield simulator to predict yields:
$$y_{i} = f(\boldsymbol{x, q; \hat{w}, \hat{s} | \alpha, \beta})$$

### Management model
We assume that a farmer has the goal of maximise revenues minus costs through optimal management:
$$\max_{\boldsymbol{x, q}} pf(\boldsymbol{x, q; \hat{w}, \hat{s} | \alpha, \beta}) - c\boldsymbol{x} - C(\boldsymbol{x} > 0)$$
- Management choices apply over a period of time $T$ for discrete time intervals (e.g. fortnightly): $t = 1, 2, ..., T$
- Planting choice: $\boldsymbol{x} = [x_{1}, x_{2}, ..., x_{T}] \quad\text{      } x_{t} &isin; \lbrace 0, 1\rbrace, \quad\text{      } \sum\limits_{t} x_{t} = 1$
- Fertiliser application choice: $\boldsymbol{q} = [q_{1}, q_{2}, ..., q_{T}] \quad\text{      } q_{t}\geq 0$
- Yield unit price: $p$
- Fertiliser application unit price: $c$
- Fertilisation event fixed cost: $C$

### Yield-managment model calibration/fitting
- Data inputs: $\boldsymbol{\hat{s}, \hat{w}}, p, c, C$
- Calibration parameters: $\boldsymbol{\alpha, \beta}$
- Outputs: $y, \boldsymbol{x, q}$

Calibration data comes from:
1. DEFRA (yields, management);
2. [Sentinel 2 athmospherically corrected (L2A)](https://registry.opendata.aws/sentinel-2/) multispectral images;
3. [Sentinel 1 granules](https://search.asf.alaska.edu/#/);
4. [ERA-5 climate data from ECMWF](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)

## Setup

### Pre-requisites

### Installation

## Usage

## Details

## Modules

## Types

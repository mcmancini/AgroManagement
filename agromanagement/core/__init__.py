# -*- coding: utf-8 -*-
# Copyright (c) 2023 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), December 2023
# ========================================================
"""
agromanagement core package initialisation file
"""
import boto3
import rasterio

from agromanagement.core.config_parser import ConfigReader
from agromanagement.utility.paths import ROOT_DIR

config_path = ROOT_DIR + "\\config.ini"
# pylint: disable=E1101
app_config = ConfigReader(config_path)
# pylint: enable=E1101

aws_session = rasterio.session.AWSSession(boto3.Session(), requester_pays=True)

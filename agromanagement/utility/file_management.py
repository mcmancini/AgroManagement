# -*- coding: utf-8 -*-
# Copyright (c) 2024 LEEP, University of Exeter (UK)
# Mattia Mancini (m.c.mancini@exeter.ac.uk), April 2024
# =====================================================
"""
File and folder management utility functions
============================================

Functions defined here:
-----------------------

unzip_all_files(input_folder, output_folder):
    find all zip archives into input_folder and unzip them
    all into output_folder

"""

import os
import zipfile


def unzip_all_files(input_folder, output_folder=None):
    """
    Find all .zip archives in 'input_folder' and extract their
    content into output_folder. If no output_folder is specified,
    the archives will be unzipped in place.

    Parameters:
    -----------
    :param input_folder (str): the location where the search for
        zip files is performed. All zip files into 'input_folder' will
        be extracted by this function
    :param output_folder (str): the location where the unzipped data will
        be stored.

    Output:
    -------
    The content extracted from the zip archives
    """
    if output_folder is None:
        output_folder = input_folder

    file_list = os.listdir(input_folder)
    zip_file_list = [
        os.path.join(input_folder, file) for file in file_list if file.endswith(".zip")
    ]

    num_files = len(zip_file_list)
    counter = 0
    for file in zip_file_list:
        with zipfile.ZipFile(file, "r") as zip_ref:
            zip_ref.extractall(output_folder)
        counter += 1
        print(f"Extracted archive {counter} of {num_files}")
    print("All archives successfully extracted!")

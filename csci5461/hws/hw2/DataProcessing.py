#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Process Rna-seq and Array data

@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: DataProcessing.py
@time: 2021/2/24 2:41 AM
"""
import click
import numpy as np
import pandas as pd


def save_data(file, survival_info, output):
    """process data by using clinic information

    :param file: filename of data
    :type file: string
    :param survival_info:  information of clinic data
    :type survival_info: pandas.DataFrame
    :param output: filename of result
    :type output: string
    """
    df = pd.read_csv(file, sep="\t", index_col=0)
    df_T = df.T
    df = df_T.join(survival_info.Group)
    df = df.loc[:, [df.columns[-1], *df.columns[:-1]]]
    df_groups = df[(df.Group == 1) | (df.Group == 2)]
    df_groups.index.name = "Sample ID"
    df_groups.reset_index().to_csv(output, index=False, sep="\t")


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("ClinicFile", type=click.Path(exists=True))
@click.argument("RNASeqFile", type=click.Path(exists=True))
@click.argument("MicroarrayFile", type=click.Path(exists=True))
def main(clinicfile, rnaseqfile, microarrayfile):
    """get parameters and manipulate files

    :param clinicfile: filename of clinic data
    :type clinicfile: string
    :param rnaseqfile: filename of Rna-seq data
    :type rnaseqfile: string
    :param microarrayfile: filename of Array data
    :type microarrayfile: string
    """

    clinic_data = pd.read_csv(clinicfile, sep="\t", index_col=["Sample ID"])
    survival_info = clinic_data.loc[
        :, ["Overall Survival Status", "Overall Survival (Months)"]
    ]

    # creat group 1
    survival_info["Group"] = np.nan
    group1_info = (survival_info["Overall Survival Status"] == "1:DECEASED") & (
        survival_info["Overall Survival (Months)"] < 36
    )
    survival_info.Group[group1_info] = 1
    # creat group 2
    group2_info = survival_info["Overall Survival (Months)"] > 36
    survival_info.Group[group2_info] = 2
    # process data
    save_data(rnaseqfile, survival_info, "SeqData.txt")
    save_data(microarrayfile, survival_info, "ArrayData.txt")

    click.echo("Processing Files Done")


if __name__ == "__main__":
    main()

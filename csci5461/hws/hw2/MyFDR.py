#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" FDR Correction

@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: MyFDR.py
@time: 2021/2/24 3:47 AM
"""
import sys

import click
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind


sys.setrecursionlimit(3000)


def ttest(datafile):
    """conduct t-test

    :param datafile: filename of data for t-test
    :type datafile: string
    :return: t-statistics, p-values and related genes
    :rtype: np.array
    """
    df = pd.read_csv(datafile, sep="\t", index_col=0)
    genes_expression = (df == 0).all(axis=0)[~(df == 0).all(axis=0)].index
    df = df[genes_expression]
    genes_id = df.columns[1:]
    group1 = df[df.Group == 1].iloc[:, 1:]
    group2 = df[df.Group == 2].iloc[:, 1:]
    t_statistics, p_values = ttest_ind(group1, group2, equal_var=False)
    return t_statistics, p_values, genes_id


def fdr_correction(pvalues):
    """fdr correction

    :param pvalues: sorted pvalues
    :type pvalues: np.array
    :return: fdr values
    :rtype: np.array
    """
    if all(sorted(pvalues) == pvalues):
        return pvalues
    else:
        temp = []
        for ind, item in enumerate(pvalues[:-1]):
            v = min(item, pvalues[ind + 1])
            temp.append(v)
        temp.append(pvalues[-1])
        return fdr_correction(np.array(temp))


def get_fdr_threshold(pvalues, numbers, genes):
    """get fdr threshold in terms of number of selected genes

    :param pvalues: unsorted and original pvalues
    :type pvalues: np.array
    :param numbers: number of genes
    :type numbers: list or tuple
    :param genes: genes that are consistent with unsorted pvalues
    :type genes: np.array
    :return: fdr threshold and genes
    :rtype: np.array
    """

    sort_index = np.argsort(pvalues)
    sort_genes = genes[sort_index]
    sort_pvalues = pvalues[sort_index]

    rank_index = np.arange(0, sort_pvalues.shape[0]) + 1

    temp_p = (sort_pvalues.shape[0] * sort_pvalues) / rank_index

    fdr = fdr_correction(temp_p)

    return fdr[np.array(numbers) - 1].tolist(), sort_genes


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("SeqData", type=click.Path(exists=True), nargs=1)
@click.argument("ArrayData", type=click.Path(exists=True), nargs=1)
@click.argument("Numbers", nargs=-1, type=click.INT)
def main(seqdata, arraydata, numbers):
    """get parameters

    :param seqdata: filename of Rna-seq data
    :type seqdata: string
    :param arraydata: filename of Array data
    :type arraydata: string
    :param numbers: numbers of selected genes
    :type numbers: list or tuple
    """
    seq_statistics, seq_pvalues, seq_genes = ttest(seqdata)
    array_statistics, array_pvalues, array_genes = ttest(arraydata)

    numbers_index = np.array(numbers) - 1

    upper_seq, _ = get_fdr_threshold(seq_pvalues, numbers_index, seq_genes)
    upper_array, _ = get_fdr_threshold(array_pvalues, numbers_index, array_genes)

    print(
        (
            f"\nSeqData:\n"
            f"Number of Genes: {numbers}\n"
            f"FDR: {upper_seq}\n"
            f"\nArrayData:\n"
            f"NUmber of Genes: {numbers}\n"
            f"FDR: {upper_array}\n"
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Bonferroni Correction

@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: MyBonferroni.py
@time: 2021/2/24 3:19 AM
"""
import click
import pandas as pd
from scipy.stats import ttest_ind


def ttest(datafile):
    """conduct t-test

    :param datafile: filename of related data
    :type datafile: string
    :return: t-statistics, p_values, and related genes
    :rtype: np.array
    """
    df = pd.read_csv(datafile, sep="\t", index_col=0)
    # remove genes with no expression
    genes_expression = (df == 0).all(axis=0)[~(df == 0).all(axis=0)].index
    df = df[genes_expression]
    genes_id = df.columns[1:]
    group1 = df[df.Group == 1].iloc[:, 1:]
    group2 = df[df.Group == 2].iloc[:, 1:]
    t_statistics, p_values = ttest_ind(group1, group2, equal_var=False)
    return t_statistics, p_values, genes_id


def bonferroni_correction(pvalues, sig, genes_id):
    """bonferroni correction for unsorted and original pvalues

    :param pvalues: original and unsorted pvalues
    :type pvalues: np.array
    :param sig: significant level
    :type sig: float
    :param genes_id: genes that are consistent with p values
    :type genes_id: np.array
    :return: selected pvalues and selected genes
    :rtype: np.array
    """
    adjusted_sig = sig / pvalues.shape[0]
    cond = pvalues < adjusted_sig
    selected_pvalue = pvalues[cond]
    selected_genes = genes_id[cond]
    return selected_pvalue, selected_genes


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("SeqData", type=click.Path(exists=True))
@click.argument("ArrayData", type=click.Path(exists=True))
@click.option("--sig", help="significant level", default=0.05, type=click.FLOAT)
def main(seqdata, arraydata, sig):
    """get parameters and conduct multi-test

    :param seqdata: filename of Rna-seq data
    :type seqdata: string
    :param arraydata: filename of Array data
    :type arraydata: string
    :param sig: significant level
    :type sig: float
    """
    # conduct t-test
    seq_statistics, seq_pvalues, seq_genes = ttest(seqdata)
    array_statistics, array_pvalues, array_genes = ttest(arraydata)
    # conduct bonferroni correction
    seq_selected_pvalues, seq_selected_genes = bonferroni_correction(
        seq_pvalues, sig, seq_genes
    )

    array_selected_pvalues, array_selected_genes = bonferroni_correction(
        array_pvalues, sig, array_genes
    )

    print(
        (
            f"\nSeqData:\n"
            f"Number of Genes: {len(seq_selected_genes)}\n"
            f"Name of Genes: {seq_selected_genes.tolist()}\n"
            f"\nArrayData:\n"
            f"NUmber of Genes: {len(array_selected_genes)}\n"
            f"Name of Genes: {array_selected_genes.tolist()}\n"
        )
    )


if __name__ == "__main__":
    main()

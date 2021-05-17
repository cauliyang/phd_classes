#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Script Used for CSCI 5461 HW4

@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: hw4.py
@time: 2021/4/18 5:54 PM
"""
import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums


def plot_degree(plot_values, name="degree_freq"):
    """ plot figure of the degree vs frequency in P1

    :param plot_values: a numpy array with two columns (degree | frequency)
    :type plot_values: np.array
    :param name: the name of figure
    :type name: string
    """
    fig, axs = plt.subplots(1, 2)

    axs = axs.flatten()

    for ind, ax in enumerate(axs):
        degree = plot_values[:, 0] if ind else np.log10(plot_values[:, 0])
        freq = plot_values[:, 1] if ind else np.log10(plot_values[:, 1])
        xlabel = "Degree" if ind else "Log Degree"
        ylbael = "Frequency" if ind else "Log Frequency"
        ax.plot(degree, freq, linewidth=2, color="r")

        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylbael, fontsize=15)

        ax.tick_params(axis="y", labelsize=12)
        ax.tick_params(axis="x", labelsize=12)

    fig.set_size_inches(14, 6)
    print(f'Saving Plot {name}.png\n')
    plt.savefig(f"{name}.png", format="png", bbox_inches="tight", dpi=300)


def get_high(number, ppi_network, df):
    """ get proteins with highest degree

    :param number: the number of proteins
    :type number: int
    :param ppi_network:  a matrix stored interactions of proteins
    :type ppi_network:   np.array
    :param df: a pandas DataFrame stored gene names and degree
    :type df: pd.DataFrame
    :return: protein ids with highest degree, interactions
    """
    high = df.sort_values(by="degree", ascending=False)[:number]
    ind_high = high.index
    gene_high = high.gene.values

    interactions = ppi_network[np.ix_(ind_high, ind_high)].sum() / 2

    return gene_high, interactions


def get_cluster_coeff(gene, ppi_network, df):
    """ calculate clustering coefficients

    :param gene: gene id
    :type gene: string
    :param ppi_network: a matrix stored interactions of proteins
    :type ppi_network: np.array
    :param df: a pandas DataFrame stored gene names and degree
    :type df: pd.DataFrame
    :return: clustering coefficient for one gene
    :rtype: float
    """
    # calculate supposed interactions
    ind = df[df.gene == gene].index.values[0]
    indx = ppi_network[ind] == 1
    supposed_interactions = (indx.sum() * (indx.sum() - 1)) / 2

    if supposed_interactions == 0:
        return 0
    # calculate actual interactions
    neibor_interactions = ppi_network[np.ix_(indx, indx)].sum() / 2

    return neibor_interactions / supposed_interactions


def plot_coeff(df, name="coeff_degree"):
    """ plot figure for clustering coefficient

    :param df: a pandas DataFrame stored gene names and degree
    :type df: pd.DataFrame
    :param name: the name of figure
    :type name: string
    """
    fig, ax = plt.subplots()

    df_temp = df.sort_values(by="coeff")

    x = df_temp["coeff"].values
    y = df_temp["degree"].values

    xlabel, ylbael = "Clustering Coefficient", "Node Degree"

    ax.plot(x, y, linewidth=2, color="navy")

    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylbael, fontsize=15)

    ax.tick_params(axis="y", labelsize=12)
    ax.tick_params(axis="x", labelsize=12)

    fig.set_size_inches(14, 6)
    print(f'Saving Plot {name}.png\n')
    plt.savefig(f"{name}.png", format="png", bbox_inches="tight", dpi=300)


def p1_solution(ppi_network, gene_names, number=5):
    """ solution for Problem 1

    :param ppi_network: a matrix stored interactions of proteins
    :type ppi_network: np.array
    :param gene_names: all gene names
    :type gene_names: list
    :param number: the number of higher proteins
    :type number: int
    :return: a pandas DataFrame stored gene names and degree
    :rtype: pd.DataFrame
    """
    # fill diagonal with 0 to remove the example in  which protein interacts with itself.
    np.fill_diagonal(ppi_network, 0)
    # get degree
    degree = ppi_network.sum(axis=0)
    # calculate frequency
    df = pd.DataFrame({"gene": gene_names, "degree": degree})
    temp_degree = df.degree.value_counts().reset_index().sort_values(by="index").to_numpy()
    temp_degree = temp_degree[temp_degree[:, 0] > 0]
    all_nodes = temp_degree[:, 1].sum()
    temp_degree[:, 1] = temp_degree[:, 1]/all_nodes
    # plot
    values_plot = temp_degree
    plot_degree(values_plot)
    # get highest proteins and interactions
    gene_high, interactions = get_high(number, ppi_network, df)

    print(
        (
            f"Proteins with the highest degrees: {gene_high}\n"
            f"Interactions: {interactions}\n"
        )
    )
    # calculate clustering coeff
    df["coeff"] = df.gene.apply(lambda x: get_cluster_coeff(x, ppi_network, df))
    # plot
    plot_coeff(df)

    return df


def p2_solution(df, gene_id: list = ['YNL110C', 'YML085C']):
    """ solution for problem 2

    :param df: a pandas DataFrame stored gene names and degree
    :type df:  pd.DataFrame
    :param gene_id: gene names for checking degree and clustering coefficient
    :type gene_id: list
    """
    result = df.set_index("gene").loc[gene_id, :]
    print(f'Node degrees and clustering coefficients for{gene_id}\n{result}\n')


def p3_solution(df):
    """ solution for Problem 3

    :param df: a pandas DataFrame stored gene names and degree
    :type df: pd.DataFrame
    :return: statistics, pvalue
    :rtype: float
    """
    group_0 = df[df.isessential == 0].degree.values
    group_1 = df[df.isessential == 1].degree.values
    s, p = ranksums(group_0, group_1)

    return s, p


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"], max_content_width=150)


@click.command(options_metavar="<options>", context_settings=CONTEXT_SETTINGS)
@click.argument("file", type=click.Path(exists=True), metavar="input")
@click.option(
    "-p1", is_flag=True, help="Run data only for problem 1", show_default=True
)
def main(file, p1):
    print('Run Program for P1\n')
    # read data
    ppi_data = np.load(file)
    ppi_network = ppi_data["ppi_network"]
    gene_names = ppi_data["gene_names"]

    # p1
    df = p1_solution(ppi_network, gene_names)

    if not p1:
        print('Run Program for P2, P3\n')
        try:
            is_essential = ppi_data['is_essential']
            essential_genes = ppi_data['essential_genes']
            df['isessential'] = is_essential
        except KeyError:
            raise SystemExit("File do not contain 'is_essential' key, you may run wrong mode! (Missing '-p1')")

        # p2
        p2_solution(df)
        # p3
        s, p = p3_solution(df)
        print(
            (
                f"Rank-Sum for node degree and protein essentiality\n"
                f"P-value: {p:e}\n"
                f"Significant Level: 0.01"
            )
        )


if __name__ == "__main__":
    main()

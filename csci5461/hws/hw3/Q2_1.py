#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: Q2_1.py
@time: 2021/3/7 5:34 PM
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.cluster import KMeans


def plot_his(labels):
    """ plot histogram of cluster size

    :param labels: labels of data
    :type labels: numpy.ndarray
    """
    k = len(set(labels))
    plt.rc('font', family='Times New Roman')
    fig, ax = plt.subplots()
    sns.countplot(y=labels + 1,
                  edgecolor=sns.color_palette("dark", 3),
                  palette='Set3',
                  ax=ax)
    ax.set_ylabel('Cluster ID', fontsize=16)
    ax.set_xlabel('Cluster Size', fontsize=16)

    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='x', labelsize=10)
    ax.set_title(f'Histogram of {k} Clusters', fontsize=18)
    fig.set_size_inches(12, 8)

    plt.savefig(f'{k}_Cluster.png', format='png', bbox_inches='tight', dpi=300)


def main(data1, k):
    """ Conduct KMeans for data1

    :param data1: data1
    :type data1: str
    :param k: k value for KMeans
    :type k: int
    """
    data1 = np.load(data1)
    selected_1000 = data1['SeqData'][:1000, ]
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(selected_1000)
    labels = kmeans.labels_
    plot_his(labels)


if __name__ == '__main__':
    pass
    # main(data1, k)

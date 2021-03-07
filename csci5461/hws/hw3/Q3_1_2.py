#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: Q3_1_2.py
@time: 2021/3/7 5:46 PM
"""
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
import numpy as np


def main(data2, k, flag):
    """ Conduct KNeighborClassifier for data2

    :param data2: data2 file
    :type data2: str
    :param k: neighbors value
    :type k: int
    :param flag: flag used to decide number of features
    :type flag: 1|0
    :return: accuracy of prediction
    :rtype: str
    """
    data2 = np.load(data2)

    if flag:
        train_data, train_label = data2['training_data'], data2['training_label']
        test_data, test_label = data2['testing_data'], data2['testing_label']
    else:
        train_data, train_label = data2['training_data'][:, :1000], data2['training_label'][:1000]
        test_data, test_label = data2['testing_data'][:, :1000], data2['testing_label'][:1000]

    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(train_data, y=train_label)
    pred = knn.predict(test_data)
    return f'{accuracy_score(test_label, pred):e}'


if __name__ == '__main__':
    pass
    # main(data2, k, flag)
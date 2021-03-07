#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: Q3_3.py
@time: 2021/3/7 5:54 PM
"""
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn.svm import LinearSVC


def main(data2):
    """ Conduct SVM for data2

    :param data2: data2 file
    :type data2: str
    :return: accuracy of prediction
    :rtype: str
    """
    data2 = np.load(data2)
    train_data, train_label = data2['training_data'][:, :1000], data2[
                                                                    'training_label'][:1000]
    test_data, test_label = data2['testing_data'][:, :1000], data2[
                                                                 'testing_label'][:1000]

    svm = LinearSVC(max_iter=5000, random_state=42)
    svm.fit(train_data, train_label)

    return f'{accuracy_score(test_label, svm.predict(test_data)):e}'


if __name__ == '__main__':
    pass
    # main(data2)

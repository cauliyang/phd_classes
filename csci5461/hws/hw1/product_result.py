#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: async_version.py
@time: 2021/1/22 3:03 PM
"""
import os
import sys

args = sys.argv[1:]


def main(args):
    if len(args) < 3:
        raise SystemExit(f"python <programe>.py sam bcf out")

    sam, bcf, out = args
    if int(sam) != 0:

        os.system(f"head -n 10 {sam} > {out}Sam.txt")
        print(f"sam file has been fetched\n")

    os.system(f"sed -n -e '30,50p' {bcf} > {out}Bcf.txt")
    print(f"bcf file has been fetched\n")


if __name__ == "__main__":

    main(args)

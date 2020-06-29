#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:01:27 2020

@author: 7ml
"""

import numpy

def say_hi(inserted):
    assert(inserted.shape==(2,2))
    print(inserted[0,0])
    released = numpy.zeros([2,3])
    print("Hello Max")
    return released
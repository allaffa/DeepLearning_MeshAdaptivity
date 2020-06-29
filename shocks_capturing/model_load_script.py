#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:36:32 2020

@author: 7ml
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json

 
def model_load():
    # load json and create model
    json_file = open('model.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into new model
    loaded_model.load_weights("model.h5")
    print("Loaded model from disk")
    return

def evaluate_model_load(input):
    # load json and create model
    json_file = open('model.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into new model
    loaded_model.load_weights("model.h5")
    print("Loaded model from disk")
    #print("input shape: ", input.shape)
    input = np.asarray(input)
    input = input.reshape(1,201)
    pred = loaded_model.predict(np.asarray(input))
    return pred

def evaluate_cnn_model_load(input):
    # load json and create model
    json_file = open('cnn_model.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into new model
    loaded_model.load_weights("cnn_model.h5")
    print("Loaded model from disk")
    print("input shape: ", input.shape)
    pred = loaded_model.predict(np.asarray(input))
    return pred






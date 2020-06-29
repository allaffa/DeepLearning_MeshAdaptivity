#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:53:49 2020

@author: 7ml
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 


#meshes = np.loadtxt('initial_meshes.txt', skiprows=0)
gradients = np.loadtxt('gradients.txt', skiprows=0)
final_deltas = np.loadtxt('final_deltas.txt', skiprows=0)
solutions = np.loadtxt('solutions.txt', skiprows=0)
adapted_meshes = np.loadtxt('adapted_meshes.txt', skiprows=0)

input_dim = solutions.shape[1]
output_dim = adapted_meshes.shape[1]
training_size = adapted_meshes.shape[0]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(solutions, adapted_meshes, test_size=0.2, random_state=42)

#build neural network
model = Sequential()
model.add(Dense(50, input_dim=input_dim, activation='relu'))
model.add(Dense(50, activation='relu'))
model.add(Dense(50, activation='relu'))
model.add(Dense(50, activation='relu'))
model.add(Dense(output_dim, activation=None))

#set-up optimizer
opt = keras.optimizers.Nadam(learning_rate=0.00001, beta_1=0.9, beta_2=0.999)

# compile the keras model
model.compile(loss='mse', optimizer=opt, metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(x_train, y_train, epochs=10000, batch_size=100, validation_split = 0.2)

#Save the trained model in a .json file
# serialize model to JSON
model_json = model.to_json()
with open("model.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model.save_weights("model.h5")
print("Saved model to disk")
 
# later...
 
# load json and create model
json_file = open('model.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights("model.h5")
print("Loaded model from disk")

results = model.evaluate(x_test, y_test, batch_size=128)
print('test loss, test acc:', results)

preds = model.predict(x_test)

num_pred = 10

uniform = np.linspace(0,1,201)

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, preds[num_pred], label = 'MLP - deep learning')
ax.plot(uniform, y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('index of the array coordinate')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, x_test[num_pred], label = 'Solution to ADR problem')
plt.xlabel('index of the array coordinate')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()
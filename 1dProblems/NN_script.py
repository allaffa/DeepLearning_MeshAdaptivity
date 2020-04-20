#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:45:37 2020

@author: 7ml (Massimiliano Lupo Pasini)
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 


initial_meshes = np.loadtxt('initial_meshes.txt', skiprows=0)
adapted_meshes = np.loadtxt('adapted_meshes.txt', skiprows=0)
#gradients = np.loadtxt('gradients.txt', skiprows=0)
solutions = np.loadtxt('solutions.txt', skiprows=0)

number_inputs = 2
input_dim = np.zeros(2)
input_dim = input_dim.astype(int)
input_dim[0] = initial_meshes.shape[1]
input_dim[1] = number_inputs
output_dim = np.zeros(2)
output_dim = input_dim.astype(int)
output_dim[0] = adapted_meshes.shape[1]
output_dim[1] = 1
data_size = solutions.shape[0]

input_data = np.zeros((data_size, input_dim[0], input_dim[1]))
adapted_meshes = np.expand_dims(adapted_meshes, axis=2)

for index in range(0,data_size):
    input_data[index, : , 0] = initial_meshes[index,:]
    input_data[index, : , 1] = solutions[index,:]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(input_data, adapted_meshes, test_size=0.2, random_state=42)

#build neural network
model = Sequential()
model.add(Dense(100, input_shape=input_dim, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(output_dim[1], activation=None))

# compile the keras model
model.compile(loss='mse', optimizer='adam', metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(x_train, y_train, epochs=50, batch_size=100, validation_split = 0.2)


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

num_pred = 97

plt.show()
uniform = np.linspace(0,1,201)
plt.plot(uniform, preds[num_pred])
plt.title('Deep learning adaptive zoning')
plt.xlabel('Coordinates in uniform mesh')
plt.ylabel('Coordinates in non-uniform mesh after adaptive zoning')
plt.savefig('DL_mesh.png')


plt.plot(uniform, y_test[num_pred])
plt.title('Standard approach for adaptive zoning')
plt.xlabel('Coordinates in uniform mesh')
plt.ylabel('Coordinates in non-uniform mesh after adaptive zoning')
plt.savefig('AdaptiveZoning_mesh.png')

plt.show()
plt.plot(x_test[num_pred, :, 0], x_test[num_pred, :, 1])
plt.title('Solution to advection-diffusion-reaction problem on uniform mesh')
plt.xlabel('Coordinates in uniform mesh')
plt.ylabel('Coordinates in non-uniform mesh after adaptive zoning')
plt.savefig('solution.png')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:45:37 2020

@author: 7ml (Massimiliano Lupo Pasini)
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense, InputLayer, Flatten
from keras.models import model_from_json
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 
from keras.regularizers import l2
from keras.optimizers import Nadam



#initial_meshes = np.loadtxt('initial_meshes.txt', skiprows=0)
#adapted_meshes = np.loadtxt('adapted_meshes.txt', skiprows=0)
#gradients = np.loadtxt('gradients.txt', skiprows=0)
#solutions = np.loadtxt('solutions.txt', skiprows=0)

initial_deltas = np.loadtxt('initial_deltas.txt', skiprows=0)
final_deltas = np.loadtxt('final_deltas.txt', skiprows=0)
gradients = np.loadtxt('gradients.txt', skiprows=0)

number_inputs = 2
input_dim = np.zeros(2)
input_dim = input_dim.astype(int)
input_dim[0] = initial_deltas.shape[1]
input_dim[1] = number_inputs
output_dim = np.zeros(2)
output_dim = input_dim.astype(int)
output_dim[0] = final_deltas.shape[1]
output_dim[1] = 1
data_size = gradients.shape[0]

input_data = np.zeros((data_size, input_dim[0], input_dim[1]))

for index in range(0,data_size):
    input_data[index, : , 0] = initial_deltas[index,:]
    input_data[index, : , 1] = gradients[index,:]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(input_data, final_deltas, test_size=0.2, random_state=42)


x_train_standard = np.zeros_like(x_train)
x_test_standard = np.zeros_like(x_test)


## We scale the data
from sklearn.preprocessing import StandardScaler

#The scaling process is based on the training set
scalers = {}
for i in range(x_train.shape[1]):
    scalers[i] = StandardScaler()
    x_train_standard[:, i, :] = scalers[i].fit_transform(x_train[:, i, :]) 

for i in range(x_test.shape[1]):
    x_test_standard[:, i, :] = scalers[i].transform(x_test[:, i, :]) 


#Set up of the neural network based on penaly and number of neurons
l2_penalty = 0.0
num_neurons = 1000

#build neural network
model = Sequential()
model.add(InputLayer((input_dim[0], input_dim[1])))
model.add(Flatten())
model.add(Dense(num_neurons, activation='relu'))
model.add(Dense(num_neurons, activation='relu'))
#model.add(Dense(num_neurons, activation='relu', kernel_regularizer=l2(l2_penalty), bias_regularizer=l2(l2_penalty)))
#model.add(Dense(num_neurons, activation='relu', kernel_regularizer=l2(l2_penalty), bias_regularizer=l2(l2_penalty)))
model.add(Dense(output_dim[0], activation='relu'))

#set-up optimizer
opt = keras.optimizers.Nadam(learning_rate=0.01, beta_1=0.9, beta_2=0.999)

# compile the keras model
model.compile(loss='mse', optimizer=opt, metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(x_train_standard, y_train, epochs=100, batch_size=100, validation_split = 0.2)


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

num_pred = 35

uniform = np.linspace(0,1,201)

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_test[num_pred, :, 0], preds[num_pred], label = 'DNN - deep learning')
ax.plot(x_test[num_pred, :, 0], y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('Coordinates in initial mesh')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()
plt.savefig('DNN_mesh.png')

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, preds[num_pred], label = 'CNN - deep learning')
ax.plot(uniform, y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('index of the array coordinate')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()


plt.show()
plt.plot(x_test[num_pred, :, 0], x_test[num_pred, :, 1])
plt.title('Solution to advection-diffusion-reaction problem')
plt.xlabel('Coordinates in initial mesh')
plt.ylabel('Value of the solution')
plt.savefig('DNN_solution.png')
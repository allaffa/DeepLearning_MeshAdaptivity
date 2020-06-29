#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 23:29:02 2020

@author: 7ml
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense, Conv1D, Flatten
from keras.models import model_from_json
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 


final_meshes = np.loadtxt('final_meshes.txt', skiprows=0)
shocks = np.loadtxt('shocks.txt', skiprows=0)
shocks = np.reshape(shocks, (shocks.shape[0], shocks.shape[1], 1))

input_dim = shocks.shape
output_dim = final_meshes.shape[1]
training_size = final_meshes.shape[0]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(shocks, final_meshes, test_size=0.2, random_state=42)


#Set up of the neural network based on penaly and number of neurons
l2_penalty = 0.0
num_neurons = 10
num_epochs = 100
batch_size = 100
activation_function = 'relu'

#set-up optimizer
opt = keras.optimizers.Nadam(learning_rate=0.0001, beta_1=0.9, beta_2=0.999)

model = Sequential()
model.add(Conv1D(filters=100, kernel_size=3, activation=activation_function, input_shape=input_dim[1:]))
model.add(Conv1D(filters=100, kernel_size=3, activation=activation_function))
#model.add(Conv1D(filters=256, kernel_size=3, activation=activation_function))
#model.add(Conv1D(filters=256, kernel_size=3, activation=activation_function))
#    model.add(Dropout(0.5))
#model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(num_neurons, activation=activation_function))
model.add(Dense(output_dim, activation=None))

# compile the keras model
model.compile(loss='mse', optimizer=opt, metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(x_train, y_train, epochs=num_epochs, batch_size=batch_size, validation_split = 0.2)

scores = model.evaluate(x_test, y_test, verbose=0)
print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))


#Save the trained model in a .json file
# serialize model to JSON
model_json = model.to_json()
with open("cnn_model.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model.save_weights("cnn_model.h5")
print("Saved model to disk")
 
# later...
 
# load json and create model
json_file = open('cnn_model.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights("cnn_model.h5")
print("Loaded model from disk")


results = model.evaluate(x_test, y_test, batch_size=128)
print('test loss, test acc:', results)

preds = model.predict(x_test)

num_pred = 35

uniform = np.linspace(0,1,201)

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, preds[num_pred], label = 'CNN - deep learning')
ax.plot(uniform, y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('index of the array coordinate')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, x_test[num_pred], label = 'Shock profile')
plt.xlabel('value of profile')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()
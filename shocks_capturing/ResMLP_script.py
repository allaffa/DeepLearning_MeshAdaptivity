#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:53:49 2020

@author: 7ml
"""

import numpy as np
import keras 
from keras.models import Sequential, Model
from keras.layers import Dense, Input, Add
from keras.models import model_from_json
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 

plt.rcParams.update({'font.size': 16})

omegas = np.loadtxt('omegas.txt', skiprows=0)
deltas = np.loadtxt('deltas.txt', skiprows=0)

input_dim = omegas.shape[1]
output_dim = deltas.shape[1]
training_size = deltas.shape[0]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(omegas, deltas, test_size=0.2, random_state=42)

#build MLP neural network
"""
model = Sequential()
model.add(Dense(100, input_dim=input_dim, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(output_dim, activation=None))
"""

input_x = Input(shape=input_dim)
x = Dense(100, activation='relu')(input_x)
x = Dense(100, activation='relu')(x)
y = Dense(100, activation='relu')(x)
addition = Add()([x, y])
x = Dense(100, activation='relu')(addition)
x = Dense(100, activation='relu')(x)
y = Dense(100, activation='relu')(x)
addition = Add()([x, y])
x = Dense(100, activation='relu')(addition)
x = Dense(100, activation='relu')(x)
y = Dense(100, activation='relu')(x)
addition = Add()([x, y])
x = Dense(100, activation='relu')(addition)
x = Dense(100, activation='relu')(x)
y = Dense(100, activation='relu')(x)
addition = Add()([x, y])
output = Dense(output_dim, activation=None)(addition)

model = Model(input_x, output)

#set-up optimizer
opt = keras.optimizers.Nadam(learning_rate=1e-4, beta_1=0.9, beta_2=0.999)

# compile the keras model
model.compile(loss='mae', optimizer=opt, metrics=['accuracy'])

# fit the keras model on the dataset
history = model.fit(x_train, y_train, epochs=5000, batch_size=100, validation_split = 0.2)

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

num_pred = 1000

uniform = np.linspace(0,1,x_train.shape[1])

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, preds[num_pred], label = 'MLP - deep learning')
ax.plot(uniform, y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('Uniform mesh node coordinates')
plt.ylabel('Adapted mesh node coordinates')
#plt.ylabel('Delta X mesh')
ax.legend()
plt.savefig('old_mesh')

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(uniform, x_test[num_pred], label = 'Shock profile')
plt.xlabel('Uniform mesh node coordinates')
plt.ylabel('Values of Shock')
#plt.ylabel('Values of Omega')
ax.legend()

"""
old_pred = model.predict(np.expand_dims(x_test[num_pred,:], axis=0))
x_test[num_pred][140:170] = 9.0
new_pred = model.predict(np.expand_dims(x_test[num_pred,:], axis=0))
"""

loss = history.history['loss']
val_loss = history.history['val_loss']

epochs = range(1,len(loss)+1)
plt.figure()
plt.plot(epochs, loss, 'b', label = 'training loss')
plt.plot(epochs, val_loss, 'r', label = 'validation loss')
plt.yscale('log')
#plt.title('Training and validation accuracy')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.draw()
plt.savefig('loss_function_plot')
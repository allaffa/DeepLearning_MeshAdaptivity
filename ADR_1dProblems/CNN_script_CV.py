#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:35:46 2020

@author: 7ml
"""

"""
Created on Tue Mar  3 11:45:37 2020

@author: 7ml (Massimiliano Lupo Pasini)
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Conv1D, Dropout, MaxPooling1D, Flatten, Dense
from keras.models import model_from_json
from keras.models import clone_model
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 
from keras.regularizers import l2
from keras.optimizers import Nadam
from sklearn.model_selection import KFold

# create a model from the weights of multiple models
def model_weight_ensemble(members, weights):
	# determine how many layers need to be averaged
	n_layers = len(members[0].get_weights())
	# create an set of average model weights
	avg_model_weights = list()
	for layer in range(n_layers):
		# collect this layer from each model
		layer_weights = np.array([model.get_weights()[layer] for model in members])
		# weighted average of weights for this layer
		avg_layer_weights = np.average(layer_weights, axis=0, weights=weights)
		# store average layer weights
		avg_model_weights.append(avg_layer_weights)
	# create a new model with the same structure
	model = clone_model(members[0])
	# set the weights in the new
	model.set_weights(avg_model_weights)
	model.compile(loss='mse', optimizer='adam', metrics=['accuracy'])
	return model



# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

# define 10-fold cross validation test harness
num_folds = 10
kfold = KFold(n_splits=num_folds, shuffle=True, random_state=seed)
cvscores = []

models = []
weights = [1/num_folds for i in range(1, num_folds+1)]

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

for index in range(0,data_size):
    input_data[index, : , 0] = initial_meshes[index,:]
    input_data[index, : , 1] = solutions[index,:]

X = input_data
Y = adapted_meshes

#Set up of the neural network based on penaly and number of neurons
l2_penalty = 0.0
num_neurons = 100
num_epochs = 100
batch_size = 100
activation_function = 'relu'

#set-up optimizer
opt = keras.optimizers.Nadam(learning_rate=0.002, beta_1=0.9, beta_2=0.999)

for train, test in kfold.split(X):
    
    model = Sequential()
    model.add(Conv1D(filters=64, kernel_size=3, activation=activation_function, input_shape=input_dim))
    model.add(Conv1D(filters=64, kernel_size=3, activation=activation_function))
#    model.add(Dropout(0.5))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(num_neurons, activation=activation_function))
    model.add(Dense(output_dim[0], activation=None))
    
    # compile the keras model
    model.compile(loss='mse', optimizer=opt, metrics=['accuracy'])

    # fit the keras model on the dataset
    model.fit(x_train, y_train, epochs=num_epochs, batch_size=batch_size, validation_split = 1.0/float(num_folds))
    
    scores = model.evaluate(X[test], Y[test], verbose=0)
    print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
    models.append(model)
    cvscores.append(scores[1] * 100)

print("%.2f%% (+/- %.2f%%)" % (np.mean(cvscores), np.std(cvscores)))

# create a new model with the weighted average of all model weights
model_averaged = model_weight_ensemble(models, weights)

# summarize the created model
model_averaged.summary()


#Save the trained model in a .json file
# serialize model to JSON
model_json = model_averaged.to_json()
with open("cnn_model.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model_averaged.save_weights("cnn_model.h5")
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


results = model_averaged.evaluate(x_test, y_test, batch_size=128)
print('test loss, test acc:', results)

preds = model_averaged.predict(x_test)

num_pred = 35

uniform = np.linspace(0,1,201)

plt.show()
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_test[num_pred, :, 0], preds[num_pred], label = 'CNN - deep learning')
ax.plot(x_test[num_pred, :, 0], y_test[num_pred], label = 'Adaptive zoning')
plt.xlabel('Coordinates in initial mesh')
plt.ylabel('Coordinates in final mesh zoning')
ax.legend()
plt.savefig('CNN_mesh.png')

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
plt.savefig('CNN_solution.png')
# iot

There are two folders with codes in them: ``Control-Based Algorithm`` Folder and ``ML-based Algorithm`` Folder

## Machine-Learning Based algorithm
* ``classifier.m`` is the implementation of ML-based rate control algorithm.
* ``NeuralNetwork.m`` is the trained model to make prediction.
* ``data_preprocessing.py`` is used to preprocess the training data.
* ``encode_model.m`` is a simple function to encode the delay model.

## Control-Based Algorithm Folder

This folder contains the Control-Based Algorithm RCA.
* ``Rate_Control_Algorithm_1.m`` contains the main process of Control-Based algorithm.
* ``compensate.m`` contains the outer loop long-term feedback codes.
* ``predict.m`` contains the SNR prediction and SNR-MCS mapping codes.

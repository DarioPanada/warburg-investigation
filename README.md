# Warburg Investigation

This repository contains code for my cancer modelling related to investigating the Warburg Effect.

## Installation

The code is designed to run under Python 2.7 . 

It is recommended to run the code inside the provided docker container. A dockerfile is provided, and the image is also available
[on dockerhub](https://hub.docker.com/r/dashma94/panaxea-warburg).

A requirements.txt file is provided specifying dependencies.

You can run the requirements.txt as `pip install -r requirements.txt`. _But_, this is [known to cause issues installing pysparse](https://github.com/usnistgov/fipy/issues/435), which is why we recommend not trying to automatically install dependencies via pip but rather installing them individually. Or, using the provided docker environment.

## Contents
* **analyzers** - Contains functions to analyze model output;
* **docker** - Contains the dockerfile with all dependencies for the code to run;
* **experiments** - Contains experiment csv files;
* **model** - Contains the model files, including all agent classes, helpers, etc.
* **requirements** - Contains experiment outputs;

## Running the Code

Your `PYTHONPATH` variable should point to the root directory of the project.

Running Main.py will run all experiments sequentially, each experiment will create a directory in `requirements` where output will be stored.

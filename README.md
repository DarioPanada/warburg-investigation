# Warburg Investigation

This repository contains code for my cancer modelling related to investigating the Warburg Effect.

## Installation

The code is designed to run under Python 2.7 . 

It is recommended to run the code inside the provided docker container. This is because the installation of all dependencies is not straightforward. A dockerfile is provided, the image is also available [on dockerhub](https://hub.docker.com/r/dashma94/panaxea-new).

A requirements.txt file is provided enumerating the dependencies. However, the installation of some modules has been failing at times so that may not work. 

You can run the setup.txt as `pip install -r requirements.txt`

## Contents

* **model** - Contains the model files, including all agent classes, helpers, etc.

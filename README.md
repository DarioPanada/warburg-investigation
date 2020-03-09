# Warburg Investigation

This repository contains code for my cancer modelling related to investigating the Warburg Effect.

## Installation

The code is designed to run under Python 2.7 . 

It is recommended to run the code inside the provided docker container. A dockerfile is provided, and the image is also available
[on dockerhub](https://hub.docker.com/r/dashma94/panaxea-warburg). The file to generate it is Dockerfile. **Note that this will *not* include the code.**

If you just want to quickly run the code, you can use Dockerfile_With_Code. This will setup a container with all dependencies and copy the code it it. The code will be copied to /warburg-investigation and the PYTHONPATH set appropriately. Then, you can just cd to the directory and run `python Main.py --pysparse`.

Dockerfile also includes some dependencies useful for cloud deployment on AWS that Dockerfile_With_Code doesn't have.

A requirements.txt file is provided specifying dependencies.

You can run the requirements.txt as `pip install -r requirements.txt`. _But_, this is [known to cause issues installing pysparse](https://github.com/usnistgov/fipy/issues/435), which is why we recommend not trying to automatically install dependencies via pip but rather installing them individually. Or, using the provided docker environment.

## Contents
* **analyzers** - Contains functions to analyze model output;
* **experiments** - Contains experiment csv files;
* **model** - Contains the model files, including all agent classes, helpers, etc.
* **reports** - Contains experiment outputs;

## Running the Code

Your `PYTHONPATH` variable should point to the root directory of the project.

Running Main.py will run all experiments sequentially, each experiment will create a directory in `reports` where output will be stored. If PySparse is available (and it is in our containers) make sure to add `--pysparse`  this *significantly* improves simulation speed.

## TODO

- Finish documenting all methods;
- Add a launcher that obtains experiments from queue and, upon completion, uploads result to a bucket;
- Add a script to automatically add mean error to each completed simulation;
- Add a "one-click" Jupyter notebook for analysis.

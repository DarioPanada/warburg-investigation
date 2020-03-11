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
* **analysis** - Contains output files generated by analyzers;
* **analyzers** - Contains functions to analyze model output;
* **aws** - Contains functions to read from and write to aws queues;
* **experiments** - Contains experiment csv files;
* **model** - Contains the model files, including all agent classes, helpers, etc.
* **reports** - Contains experiment outputs;

## Running the Code

Your `PYTHONPATH` variable should point to the root directory of the project.

Running Main.py will run all experiments sequentially, each experiment will create a directory in `reports` where output will be stored. If PySparse is available (and it is in our containers) make sure to add `--pysparse`  this *significantly* improves simulation speed.

## Running the Experiments End-to-End

This is a quick guide on how to run the experiments end-to-end.

### Simulation Setup

The experiments file is in ./experiments/experiments_warburg.csv . For each simulation, this details all parameter values.

### Running the Simulations

#### Running Locally

This is probably only useful for demo purposes as running all experiments on a personal computer would take a lot of time.

Running `python Main.py` (or `python Main.py --pysparse` where pysparse is available) will run
all of the experiments sequentially. Outputs are saved in `./reports`, where each experiment will create a directory named as itself. The directories will contain an `./imgs` subdirectory which will include some auto-generated graphs detailing the progression.

Perhaps more importantly, the report directory will contain a *.pickle* file. This is a stripped down serialization of the model's state at simulation ends. A lot of the features (Eg: Schedule, Environments, etc.) may have been removed to save memory. **But**, it will include the output key which allows to retrieve properties such as number of agents at each epoch, hif distributions, and essentially anything else that has been stored by helpers.

#### Running on AWS

TODO

### Running the Analysis

There are two useful functions in `./analyzers/ModelErrorFunctions.py`. 

`add_ame_to_experiments` will create a csv file identical to the input experiment file but with an added column for average mean error (ame). Only experiments which have a matching report are included.

`add_error_series` will create a csv file identical to the input experiment file but with an added column for each epoch the simulation was intended to run for, whose value will be the error at such epoch. Only experiments which have a matching report are included.

A stub is provided, so the script should be runnable with minimal effort.

A jupyter notebook, `./analyzers/Warburg Analysis.ipynb` is provided. In order to run, **this requires csv files generated by the aforementioned scripts to be available**. Running it end-to-end, it allows to visualize the impact of pWarburgSwitch and minimMumOxygenConcentration on average mean error, as well as on the average direction of error.

## TODO

- Finish documenting all methods;
- Add a launcher that obtains experiments from queue and, upon completion, uploads result to a bucket;

FROM continuumio/anaconda2:5.3.0
USER root
RUN apt-get install git
ENV PYTHONPATH /app #equivalent to saying export PYTHONPATH=/app
ENV AWS_DEFAULT_REGION us-east-2

# setting up venv

# Installing necessary packages
RUN conda install -c conda-forge pympler
RUN pip install boto3 
RUN pip install awscli
RUN git clone https://github.com/usnistgov/fipy.git
RUN cd fipy && python setup.py install
RUN conda install -c conda-forge pyhamcrest
RUN conda install -c anaconda requests
RUN pip install simplejson
RUN pip install argparse 
RUN pip install setuptools
RUN apt-get update
RUN apt-get install build-essential -y
RUN pip install -e git+http://git.code.sf.net/p/pysparse/git#egg=PySparse

WORKDIR /app
COPY . /app


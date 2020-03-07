FROM continuumio/anaconda2:5.3.0
USER root
RUN apt-get install git
ENV PYTHONPATH /app #equivalent to saying export PYTHONPATH=/app
ENV AWS_DEFAULT_REGION us-east-2

# setting up venv

RUN pip install --upgrade pip

# Installing necessary packages
RUN apt-get update
RUN apt-get install build-essential -y
RUN pip install -e git+http://git.code.sf.net/p/pysparse/git#egg=PySparse
RUN pip install panaxea
RUN pip install pympler
RUN pip install boto3 
RUN pip install awscli
RUN pip install future
RUN pip install fipy
RUN pip install pyhamcrest
RUN pip install requests
RUN pip install simplejson
RUN pip install argparse 
RUN pip install setuptools

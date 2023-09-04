# syntax=docker/dockerfile:1
ARG PY_VERSION=3.9
FROM python:${PY_VERSION} as base

# Set the working directory inside the container
WORKDIR /gene-planet

COPY . .

# Copy the requirements file to the container
COPY requirements.txt .

# Install the app dependencies
RUN pip install --no-cache-dir -r requirements.txt

#ENTRYPOINT ["python", "runner.py"]
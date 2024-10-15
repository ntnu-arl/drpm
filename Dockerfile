# Use an official Ubuntu base image
FROM ubuntu:20.04

# Set the working directory in the container
WORKDIR /example

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the local project files to the container's workspace
ADD . /example

# Build the C++ project
RUN mkdir /example/build && cd /example/build \
    && cmake .. \
    && make

# Command to run when starting the container
CMD ["./build/example"]
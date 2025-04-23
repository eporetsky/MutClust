FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.9 \
    python3.9-dev \
    python3-pip \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set up Python 3.9 as default
RUN ln -sf /usr/bin/python3.9 /usr/bin/python3 && \
    ln -sf /usr/bin/python3.9 /usr/bin/python && \
    ln -sf /usr/bin/pip3 /usr/bin/pip

# Install pytest
RUN pip install pytest

# Create a working directory
WORKDIR /app

# Copy the project files
COPY . /app/

# Install MutClust and its dependencies
RUN pip install --upgrade pip && \
    pip install .

# Create a directory for mounting data
RUN mkdir /data

# Set the entrypoint
ENTRYPOINT ["mutclust"] 
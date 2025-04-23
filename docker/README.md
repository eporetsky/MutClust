# MutClust Docker Container

This directory contains the Docker setup for running MutClust in a containerized environment.

## Building the Container

To build the Docker container, run the following command from the root directory of the MutClust repository:

```bash
docker build -t mutclust .
```

## Running MutClust

To run MutClust using the Docker container, you can use the following command:

```bash
docker run -v /path/to/your/data:/data mutclust [mutclust arguments]
```

Replace `/path/to/your/data` with the path to your data directory on your host machine. The container will mount this directory to `/data` inside the container.

## Example Usage

```bash
# Example: Run MutClust with your data mounted at /data
docker run -v /path/to/your/data:/data mutclust --input /data/your_input_file.csv --output /data/results
```

## Notes

- The container uses Ubuntu 18.04 as the base image
- Python 3.9 is installed and configured as the default Python version
- All MutClust dependencies are automatically installed during the build process
- The container's entrypoint is set to the `mutclust` command, so you can directly pass MutClust arguments to the container 
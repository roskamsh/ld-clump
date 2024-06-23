# Use the official Ubuntu image as a base image
FROM ubuntu:latest

# Install wget and gunzip
RUN apt-get update && apt-get install -y wget gzip && apt-get clean

# Set the entrypoint to an interactive shell
ENTRYPOINT ["/bin/bash"]
# Docker file for creating continer for plink 1.9
# Use a base image with a compatible Linux distribution
FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && \
    apt-get install -y \
        wget \
        unzip \
        && rm -rf /var/lib/apt/lists/*

# Download PLINK
WORKDIR /tmp
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink && \
    rm plink_linux_x86_64_20231211.zip

# Add PLINK directory to PATH
ENV PATH="/usr/local/bin:${PATH}"
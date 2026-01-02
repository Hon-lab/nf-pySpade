# Dockerfile for nf-pySpade pipeline
# This container includes pySpade and the custom helper scripts for the Nextflow pipeline

FROM igvf/pyspade:pyspade_0.1.7

# Metadata
LABEL maintainer="nf-pySpade"
LABEL description="pySpade container with custom Nextflow pipeline scripts"
LABEL version="0.1.7-nf"

# Install any additional dependencies if needed
# RUN pip install --no-cache-dir <additional-packages>

# Copy the custom Python helper scripts into the container
COPY script/ /opt/nf-pyspade/script/

# Make scripts executable
RUN chmod +x /opt/nf-pyspade/script/*.py

# Add scripts to PATH
ENV PATH="/opt/nf-pyspade/script:${PATH}"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]

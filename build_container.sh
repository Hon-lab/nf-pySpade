#!/bin/bash
# Build script for nf-pySpade container images
# This script builds both Docker and Singularity containers with the custom scripts

set -e

# Version information
VERSION="0.1.7-nf"
DOCKER_IMAGE="nf-pyspade"
SINGULARITY_IMAGE="nf-pyspade_${VERSION}.sif"

echo "=========================================="
echo "nf-pySpade Container Build Script"
echo "Version: ${VERSION}"
echo "=========================================="
echo ""

# Function to build Docker image
build_docker() {
    echo "Building Docker image..."
    docker build -t ${DOCKER_IMAGE}:${VERSION} -t ${DOCKER_IMAGE}:latest .
    
    if [ $? -eq 0 ]; then
        echo "✓ Docker image built successfully: ${DOCKER_IMAGE}:${VERSION}"
        echo ""
        echo "To push to Docker Hub:"
        echo "  docker tag ${DOCKER_IMAGE}:${VERSION} <your-username>/${DOCKER_IMAGE}:${VERSION}"
        echo "  docker push <your-username>/${DOCKER_IMAGE}:${VERSION}"
    else
        echo "✗ Docker build failed"
        return 1
    fi
}

# Function to build Singularity image
build_singularity() {
    echo "Building Singularity image..."
    
    if ! command -v singularity &> /dev/null; then
        echo "⚠ Singularity not found. Skipping Singularity build."
        echo "  Install Singularity or build on a system with Singularity available."
        return 1
    fi
    
    singularity build ${SINGULARITY_IMAGE} Singularity.def
    
    if [ $? -eq 0 ]; then
        echo "✓ Singularity image built successfully: ${SINGULARITY_IMAGE}"
        echo ""
        echo "To use the image:"
        echo "  singularity exec ${SINGULARITY_IMAGE} <command>"
    else
        echo "✗ Singularity build failed"
        return 1
    fi
}

# Function to build from Docker to Singularity
build_singularity_from_docker() {
    echo "Building Singularity image from Docker..."
    
    if ! command -v singularity &> /dev/null; then
        echo "⚠ Singularity not found."
        return 1
    fi
    
    # Build from local Docker image
    singularity build ${SINGULARITY_IMAGE} docker-daemon://${DOCKER_IMAGE}:${VERSION}
    
    if [ $? -eq 0 ]; then
        echo "✓ Singularity image built from Docker successfully: ${SINGULARITY_IMAGE}"
    else
        echo "✗ Singularity build from Docker failed"
        return 1
    fi
}

# Main script
case "${1}" in
    docker)
        build_docker
        ;;
    singularity)
        build_singularity
        ;;
    both)
        build_docker
        echo ""
        build_singularity_from_docker
        ;;
    *)
        echo "Usage: $0 {docker|singularity|both}"
        echo ""
        echo "Options:"
        echo "  docker       - Build Docker image only"
        echo "  singularity  - Build Singularity image from definition file"
        echo "  both         - Build Docker image, then convert to Singularity"
        echo ""
        echo "Examples:"
        echo "  $0 docker"
        echo "  $0 singularity"
        echo "  $0 both"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Build process complete!"
echo "=========================================="

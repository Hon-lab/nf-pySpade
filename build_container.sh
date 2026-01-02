#!/bin/bash
# Build script for nf-pySpade container images
# This script builds both Docker/Podman and Singularity containers with the custom scripts
# Supports rootless builds with Podman

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

# Detect container runtime (prefer podman for rootless, fallback to docker)
detect_runtime() {
    if command -v podman &> /dev/null; then
        CONTAINER_RUNTIME="podman"
        echo "Using Podman (rootless compatible)"
    elif command -v docker &> /dev/null; then
        CONTAINER_RUNTIME="docker"
        echo "Using Docker"
    else
        echo "✗ Neither Podman nor Docker found. Please install one."
        exit 1
    fi
    echo ""
}

# Function to build Docker/Podman image
build_docker() {
    detect_runtime
    
    echo "Building container image with ${CONTAINER_RUNTIME}..."
    ${CONTAINER_RUNTIME} build -t ${DOCKER_IMAGE}:${VERSION} -t ${DOCKER_IMAGE}:latest .
    
    if [ $? -eq 0 ]; then
        echo "✓ Container image built successfully: ${DOCKER_IMAGE}:${VERSION}"
        echo ""
        echo "To push to registry:"
        if [ "$CONTAINER_RUNTIME" = "podman" ]; then
            echo "  ${CONTAINER_RUNTIME} tag ${DOCKER_IMAGE}:${VERSION} <registry>/<username>/${DOCKER_IMAGE}:${VERSION}"
            echo "  ${CONTAINER_RUNTIME} push <registry>/<username>/${DOCKER_IMAGE}:${VERSION}"
            echo ""
            echo "Example for Docker Hub:"
            echo "  ${CONTAINER_RUNTIME} tag ${DOCKER_IMAGE}:${VERSION} docker.io/<username>/${DOCKER_IMAGE}:${VERSION}"
            echo "  ${CONTAINER_RUNTIME} push docker.io/<username>/${DOCKER_IMAGE}:${VERSION}"
            echo ""
            echo "Example for GitHub Container Registry:"
            echo "  ${CONTAINER_RUNTIME} tag ${DOCKER_IMAGE}:${VERSION} ghcr.io/<username>/${DOCKER_IMAGE}:${VERSION}"
            echo "  ${CONTAINER_RUNTIME} push ghcr.io/<username>/${DOCKER_IMAGE}:${VERSION}"
        else
            echo "  docker tag ${DOCKER_IMAGE}:${VERSION} <your-username>/${DOCKER_IMAGE}:${VERSION}"
            echo "  docker push <your-username>/${DOCKER_IMAGE}:${VERSION}"
        fi
    else
        echo "✗ Container build failed"
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

# Function to build from Docker/Podman to Singularity
build_singularity_from_docker() {
    detect_runtime
    echo "Building Singularity image from ${CONTAINER_RUNTIME}..."
    
    if ! command -v singularity &> /dev/null; then
        echo "⚠ Singularity not found."
        return 1
    fi
    
    # Build from local container image
    if [ "$CONTAINER_RUNTIME" = "podman" ]; then
        # Podman uses podman-daemon protocol
        singularity build ${SINGULARITY_IMAGE} docker-daemon://${DOCKER_IMAGE}:${VERSION}
    else
        singularity build ${SINGULARITY_IMAGE} docker-daemon://${DOCKER_IMAGE}:${VERSION}
    fi
    
    if [ $? -eq 0 ]; then
        echo "✓ Singularity image built from ${CONTAINER_RUNTIME} successfully: ${SINGULARITY_IMAGE}"
    else
        echo "✗ Singularity build from ${CONTAINER_RUNTIME} failed"
        return 1
    fi
}

# Function to build Podman image specifically (for explicit podman usage)
build_podman() {
    CONTAINER_RUNTIME="podman"
    
    if ! command -v podman &> /dev/null; then
        echo "✗ Podman not found. Please install Podman."
        echo "  Installation: https://podman.io/getting-started/installation"
        return 1
    fi
    
    echo "Building container image with Podman (rootless)..."
    podman build -t ${DOCKER_IMAGE}:${VERSION} -t ${DOCKER_IMAGE}:latest .
    
    if [ $? -eq 0 ]; then
        echo "✓ Podman image built successfully: ${DOCKER_IMAGE}:${VERSION}"
        echo ""
        echo "To push to registry:"
        echo "  podman tag ${DOCKER_IMAGE}:${VERSION} <registry>/<username>/${DOCKER_IMAGE}:${VERSION}"
        echo "  podman push <registry>/<username>/${DOCKER_IMAGE}:${VERSION}"
        echo ""
        echo "Example for Docker Hub:"
        echo "  podman tag ${DOCKER_IMAGE}:${VERSION} docker.io/<username>/${DOCKER_IMAGE}:${VERSION}"
        echo "  podman push docker.io/<username>/${DOCKER_IMAGE}:${VERSION}"
        echo ""
        echo "Example for GitHub Container Registry (GHCR):"
        echo "  podman tag ${DOCKER_IMAGE}:${VERSION} ghcr.io/<username>/${DOCKER_IMAGE}:${VERSION}"
        echo "  podman push ghcr.io/<username>/${DOCKER_IMAGE}:${VERSION}"
    else
        echo "✗ Podman build failed"
        return 1
    fi
}

# Main script
case "${1}" in
    docker)
        build_docker
        ;;
    podman)
        build_podman
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
        echo "Usage: $0 {docker|podman|singularity|both}"
        echo ""
        echo "Options:"
        echo "  docker       - Build container image (auto-detects podman/docker)"
        echo "  podman       - Build container image explicitly with Podman (rootless)"
        echo "  singularity  - Build Singularity image from definition file"
        echo "  both         - Build container image, then convert to Singularity"
        echo ""
        echo "Examples:"
        echo "  $0 docker      # Auto-detects and uses podman if available"
        echo "  $0 podman      # Explicitly use podman for rootless build"
        echo "  $0 singularity"
        echo "  $0 both"
        echo ""
        echo "Note: Podman is recommended for rootless (non-root) container builds."
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Build process complete!"
echo "=========================================="

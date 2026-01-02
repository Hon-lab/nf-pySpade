# Container Build Guide for nf-pySpade

This directory contains files for building custom container images that include the nf-pySpade helper scripts along with the base pySpade package.

## Files

- **Dockerfile** - Docker container definition
- **Singularity.def** - Singularity container definition
- **build_container.sh** - Automated build script

## Quick Start

### Option 1: Using the Build Script

```bash
# Build Docker image
./build_container.sh docker

# Build Singularity image
./build_container.sh singularity

# Build both (Docker first, then convert to Singularity)
./build_container.sh both
```

### Option 2: Manual Docker Build

```bash
# Build the Docker image
docker build -t nf-pyspade:0.1.7-nf .

# Test the image
docker run --rm nf-pyspade:0.1.7-nf pySpade --help

# Tag for Docker Hub (optional)
docker tag nf-pyspade:0.1.7-nf <your-username>/nf-pyspade:0.1.7-nf

# Push to Docker Hub (optional)
docker push <your-username>/nf-pyspade:0.1.7-nf
```

### Option 3: Manual Singularity Build

```bash
# Build from definition file
singularity build nf-pyspade_0.1.7-nf.sif Singularity.def

# OR build from Docker image
singularity build nf-pyspade_0.1.7-nf.sif docker://igvf/pyspade:pyspade_0.1.7

# Test the image
singularity exec nf-pyspade_0.1.7-nf.sif pySpade --help
```

## Container Contents

The custom container includes:

1. **Base Image**: `igvf/pyspade:pyspade_0.1.7`
   - pySpade package for differential expression analysis
   - All required Python dependencies
   - GPU support (CUDA-enabled)

2. **Custom Helper Scripts** (in `/opt/nf-pyspade/script/`):
   - `calculate_FDR.py` - Calculate FDR and significance score cutoffs
   - `filtered_local_df.py` - Filter local differential expression results
   - `find_DErand_range.py` - Determine cell number distribution for DErand
   - `randomized_sgrna.py` - Generate randomized sgRNA matrix for FDR estimation

3. **Environment**:
   - Scripts added to PATH for easy execution
   - All scripts marked as executable

## Using the Container with Nextflow

### Docker

Update `nextflow.config`:
```groovy
process.container = 'docker://<your-username>/nf-pyspade:0.1.7-nf'
docker.enabled = true
```

Or keep Singularity enabled and it will auto-convert:
```groovy
process.container = 'docker://<your-username>/nf-pyspade:0.1.7-nf'
singularity.enabled = true
```

### Singularity (Local File)

If you built a local Singularity image:

```groovy
process.container = './nf-pyspade_0.1.7-nf.sif'
singularity.enabled = true
```

### Singularity (Docker Hub)

Singularity can pull directly from Docker Hub:

```groovy
process.container = 'docker://<your-username>/nf-pyspade:0.1.7-nf'
singularity.enabled = true
```

## Updating the Pipeline to Use Custom Container

If using a custom container with scripts included, you may want to update the script paths in `main.nf`:

**Current** (scripts expected in output directory):
```bash
$outdir/script/randomized_sgrna.py -s $outdir/Singlet_sgRNA_df.h5 -o ...
```

**Updated** (scripts in container PATH):
```bash
randomized_sgrna.py -s $outdir/Singlet_sgRNA_df.h5 -o ...
```

Or reference by full path:
```bash
/opt/nf-pyspade/script/randomized_sgrna.py -s $outdir/Singlet_sgRNA_df.h5 -o ...
```

## Building on HPC Systems

Many HPC systems don't allow Docker but support Singularity:

### Method 1: Build locally, transfer to HPC
```bash
# On local machine with Docker
./build_container.sh docker
singularity build nf-pyspade_0.1.7-nf.sif docker-daemon://nf-pyspade:0.1.7-nf

# Transfer to HPC
scp nf-pyspade_0.1.7-nf.sif username@hpc-system:/path/to/containers/
```

### Method 2: Build on HPC with Singularity
```bash
# On HPC system with Singularity
singularity build nf-pyspade_0.1.7-nf.sif Singularity.def
```

### Method 3: Pull from Docker Hub on HPC
```bash
# On HPC system
singularity pull docker://<your-username>/nf-pyspade:0.1.7-nf
```

## Testing the Container

### Test pySpade Commands
```bash
# Docker
docker run --rm nf-pyspade:0.1.7-nf pySpade --help

# Singularity
singularity exec nf-pyspade_0.1.7-nf.sif pySpade --help
```

### Test Helper Scripts
```bash
# Docker
docker run --rm nf-pyspade:0.1.7-nf calculate_FDR.py --help

# Singularity
singularity exec nf-pyspade_0.1.7-nf.sif calculate_FDR.py --help
```

### Test GPU Support (if available)
```bash
# Docker
docker run --gpus all --rm nf-pyspade:0.1.7-nf python -c "import torch; print(torch.cuda.is_available())"

# Singularity
singularity exec --nv nf-pyspade_0.1.7-nf.sif python -c "import torch; print(torch.cuda.is_available())"
```

## Troubleshooting

### Docker Build Issues

**Problem**: Cannot find script files
```
COPY failed: file not found in build context
```
**Solution**: Ensure you're running `docker build` from the repository root where the `script/` directory exists.

**Problem**: Permission denied
```
docker: Got permission denied while trying to connect to the Docker daemon socket
```
**Solution**: Add your user to the docker group: `sudo usermod -aG docker $USER` (then log out and back in)

### Singularity Build Issues

**Problem**: Singularity not found
```
singularity: command not found
```
**Solution**: Build on a system with Singularity installed, or use Docker and transfer the image.

**Problem**: Build requires root/sudo
**Solution**: Use `singularity build --fakeroot` if available, or build on a system where you have appropriate permissions.

### Container Usage Issues

**Problem**: Scripts not found in PATH
**Solution**: Use full path `/opt/nf-pyspade/script/<script>.py` or check that ENV PATH is set correctly in the container.

**Problem**: Permission denied for scripts
**Solution**: Rebuild container ensuring `chmod +x` is applied in the Dockerfile/Singularity.def.

## Version History

- **0.1.7-nf** - Initial custom container with nf-pySpade helper scripts
  - Based on igvf/pyspade:pyspade_0.1.7
  - Includes GPU optimization support
  - Added custom helper scripts for Nextflow pipeline

## Support

For issues related to:
- **pySpade package**: See https://github.com/Hon-lab/pySpade
- **nf-pySpade pipeline**: See repository issues
- **Container build**: Check this guide and the build script comments

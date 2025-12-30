# GPU Acceleration Implementation Notes

## Important Considerations

### 1. PySpade GPU Support

The implementation adds the `--use-gpu` flag to the following pySpade commands:
- `pySpade DEobs --use-gpu`
- `pySpade DErand --use-gpu`

**Critical:** This assumes that pySpade version 0.1.7 (or the version in the container) supports the `--use-gpu` flag for GPU-accelerated hypergeometric tests.

#### If the flag is NOT supported:

You have two options:

**Option A: Remove the GPU flag (use GPU for parallelism only)**
Remove the `--use-gpu` flag from the three processes:
- pySpadeDEobs (line ~211 in main.nf)
- pySpadeDEobsFDR (line ~241 in main.nf)  
- pySpadeDErand (line ~277 in main.nf)

The GPU resources will still be allocated and can be used if pySpade internally detects CUDA, but the explicit flag won't be passed.

**Option B: Wait for pySpade container update**
Contact the pySpade developers to add GPU support for hypergeometric tests in a future version.

### 2. GPU Resource Allocation

The processes request GPU resources using:
```groovy
queue 'GPU,GPUp40,GPUp100,256GB,256GBv1,384GB,512GB'
clusterOptions = '--gres=gpu:1'
```

This means:
- Jobs will preferentially try the GPU partitions first (GPU, GPUp40, GPUp100)
- If GPU partitions are unavailable or full, they'll fall back to memory-based queues
- Each task requests 1 GPU (`--gres=gpu:1`)

**To verify GPU partitions exist:**
```bash
sinfo -p GPU,GPUp40,GPUp100
```

If GPU partitions don't exist, remove them from the queue list.

### 3. Container Compatibility

The container is now pulled from Docker Hub:
```
docker://igvf/pyspade:pyspade_0.1.7
```

**Verify the container includes:**
- GPU support (CUDA libraries if needed)
- The pySpade commands: DEobs, DErand, process, fc, local, global, manhattan
- Python scripts are compatible with the container environment

**To test locally:**
```bash
singularity pull docker://igvf/pyspade:pyspade_0.1.7
singularity exec pyspade_pyspade_0.1.7.sif pySpade --help
```

### 4. Script Paths

The pipeline uses `$outdir/script/` for Python helper scripts:
- randomized_sgrna.py
- find_DErand_range.py
- calculate_FDR.py
- filtered_local_df.py

**These scripts must be:**
- Present in the output directory OR
- Mounted/accessible inside the container OR
- Pre-installed in the container

Current implementation assumes they're in `$outdir/script/` which gets mounted via Singularity's autoMounts.

### 5. Working Directory Usage

All processes now include:
```bash
cd ${workflow.workDir}
```

This changes to the Nextflow work directory before execution. **However**, the actual computation may not need this since:
- Input files are staged by Nextflow
- Output paths are absolute (`$outdir/...`)

If processes fail with "cannot find file" errors, the `cd ${workflow.workDir}` line may need to be removed.

## Testing Recommendations

### 1. Test without GPU first
Comment out or remove:
- `clusterOptions = '--gres=gpu:1'`
- `--use-gpu` flags
- GPU partitions from queue lists

Run a small test to verify the pipeline works with the new container and configuration.

### 2. Test GPU allocation
On a GPU node, verify CUDA is available:
```bash
srun -p GPU --gres=gpu:1 --pty bash
# or use GPUp40 or GPUp100 partition
srun -p GPUp40 --gres=gpu:1 --pty bash
nvidia-smi
singularity exec docker://igvf/pyspade:pyspade_0.1.7 python -c "import torch; print(torch.cuda.is_available())"
```

### 3. Test resume functionality
1. Run the pipeline
2. Cancel it mid-execution (Ctrl+C)
3. Run again with `-resume`
4. Verify it continues from the last completed task

### 4. Verify reports
After a complete run, check:
- `report.html` - opens in browser, shows all tasks
- `timeline.html` - shows execution timeline
- `trace.txt` - contains detailed metrics
- `dag.svg` - shows workflow structure

## Known Limitations

1. **GPU flag assumption**: The `--use-gpu` flag may not exist in pySpade 0.1.7
2. **CUDA compatibility**: Container must include CUDA libraries matching the GPU nodes
3. **Singularity binding**: Scripts in `$outdir/script/` must be accessible to container
4. **GPU queue**: The 'GPU' partition must exist in SLURM configuration

## Rollback Instructions

To revert to the original configuration:

1. **Container**: Change back to local .sif file:
   ```groovy
   process.container = './pyspade_v0150.sif'
   ```

2. **Remove GPU requests**: Delete these lines from GPU processes:
   ```groovy
   clusterOptions = '--gres=gpu:1'
   # and remove 'GPU' from queue lists
   # and remove --use-gpu flags
   ```

3. **Remove work directory changes**: Delete:
   ```bash
   cd ${workflow.workDir}
   ```

4. **Disable reports** (optional): Set to false in nextflow.config:
   ```groovy
   report.enabled = false
   timeline.enabled = false
   trace.enabled = false
   ```

5. **Disable resume** (optional): Remove from nextflow.config:
   ```groovy
   resume = true
   ```

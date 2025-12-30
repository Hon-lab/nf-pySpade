# nf-pySpade Pipeline Optimization Notes

## Changes Made

### 1. GPU Acceleration for Hypergeometric Tests
- **DEobs and DErand processes** now request GPU resources from SLURM
- Added `clusterOptions = '--gres=gpu:1'` to request one GPU per task
- Added GPU queue to process queue list: `'GPU,256GB,256GBv1,384GB,512GB'`
- Added `--use-gpu` flag to pySpade DEobs and DErand commands
- Set `CUDA_VISIBLE_DEVICES` environment variable in GPU processes

**Note:** GPU acceleration requires that the pySpade package (inside the container) supports the `--use-gpu` flag. If this flag is not supported in the current version, the processes will fall back to CPU execution or may need the flag removed.

### 2. Crash Recovery and Resume Capability
- **Enabled `-resume` flag** in the pipeline execution script
- Set `resume = true` in nextflow.config
- Added `publishDir` directives to all processes to preserve outputs
- Configured `mode: 'copy', overwrite: false` to avoid re-executing completed tasks
- All intermediate files are cached in the `work/` directory

**To resume from a crash:**
```bash
nextflow run main.nf -resume
```

### 3. Working Directory Configuration
- Set `workDir = './work'` in nextflow.config
- All processes change to `${workflow.workDir}` before executing commands
- This ensures all temporary files and caching happen in the `work/` folder

### 4. Execution Reports and Metrics
Enabled comprehensive reporting to track time and compute resources:
- **report.html** - Execution report with task statistics
- **timeline.html** - Visual timeline of task execution
- **trace.txt** - Detailed metrics for each task including:
  - CPU usage (%), memory usage (%), RSS, VMem
  - Duration, realtime, queue time
  - Submit, start, complete timestamps
  - Read/write bytes and syscalls
- **dag.svg** - Workflow DAG visualization

### 5. Container Updates
- **Updated to latest container:** `docker://igvf/pyspade:pyspade_0.1.7`
- Container is automatically pulled from Docker Hub if not found locally
- Set `singularity.pullTimeout = '60 min'` for large container pulls
- Added `singularity.autoMounts = true` for automatic directory mounting

## Usage

### Running the Pipeline

1. **First run or after changes:**
   ```bash
   sbatch log.run_nextflow.sh
   ```

2. **Resume after crash:**
   The pipeline automatically resumes from the last completed step due to the `-resume` flag in the script.

3. **Manual execution with resume:**
   ```bash
   nextflow run main.nf -resume
   ```

### Viewing Reports

After execution completes:
- Open `report.html` in a browser for overall statistics
- Open `timeline.html` for execution timeline
- Check `trace.txt` for detailed per-task metrics
- View `dag.svg` for workflow structure

### GPU Availability

The pipeline will try to use GPU resources when available. If GPUs are not available in the requested queues, SLURM will queue the jobs until GPU resources become available, or they will execute on non-GPU nodes if the GPU request cannot be fulfilled.

To run without GPU acceleration (if needed), modify the processes to remove:
- The GPU queue from the queue list
- The `clusterOptions = '--gres=gpu:1'` line
- The `--use-gpu` flag from pySpade commands

## File Organization

```
.
├── work/                    # Nextflow working directory (cached tasks)
├── DEobs/                   # DEobs output (published from work/)
├── DErand/                  # DErand output (published from work/)
├── FDR/                     # FDR analysis output
├── Manhattan_plots/         # Final plots and filtered results
├── report.html              # Execution report
├── timeline.html            # Timeline visualization
├── trace.txt               # Detailed metrics
└── dag.svg                 # Workflow DAG
```

## Troubleshooting

### Container Pull Issues
If the container fails to pull:
1. Check internet connectivity from compute nodes
2. Verify Docker Hub accessibility
3. Manually pull using: `singularity pull docker://igvf/pyspade:pyspade_0.1.7`

### GPU Issues
If GPU jobs fail:
1. Verify GPU queue availability: `sinfo -p GPU`
2. Check GPU allocation: `squeue -p GPU`
3. Test GPU access in interactive session: `salloc -p GPU --gres=gpu:1`

### Resume Not Working
If `-resume` doesn't work as expected:
1. Don't delete the `work/` directory - it contains cached results
2. Check that `.nextflow/` directory exists
3. Ensure file timestamps haven't changed
4. Review `.nextflow.log` for errors

### Performance Monitoring
Check resource usage in `trace.txt`:
- High `%cpu` but low `realtime` - CPU bound, good performance
- Low `%cpu` - I/O bound or waiting
- High `peak_rss` - memory intensive tasks
- Compare `duration` vs `realtime` to see queuing delays

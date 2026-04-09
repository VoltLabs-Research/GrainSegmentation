# GrainSegmentation

`GrainSegmentation` segments grains from the current frame using PTM-backed local structure information.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- GrainSegmentation
```

## CLI

Usage:

```bash
grain-segmentation <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--rmsd <float>` | No | RMSD threshold for PTM. | `0.1` |
| `--minGrainAtomCount <int>` | No | Minimum atoms per grain. | `100` |
| `--adoptOrphanAtoms <true\|false>` | No | Adopt orphan atoms into neighboring grains. | `true` |
| `--handleCoherentInterfaces <true\|false>` | No | Handle coherent interfaces specially. | `true` |
| `--outputBonds` | No | Export neighbor bonds. | `false` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |

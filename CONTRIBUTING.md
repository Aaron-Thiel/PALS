# Contributing to PALS

Contributions to PALS are welcome. This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository
2. Clone your fork locally
3. Create a new branch for your feature or fix
4. Make your changes
5. Test your changes
6. Submit a pull request

## Development Setup

### Prerequisites

- Nextflow >= 21.04.0
- Docker
- Conda (for visualization processes)

### Setting up the development environment

```bash
git clone <your-fork-url>
cd nextflow

# Download required databases (see info/DATABASE_SETUP.md)
./info/download_databases.sh

# Run the test profile to verify everything works
nextflow run assembly.nf -c assembly.config -profile test
```

## Code Style

### Nextflow conventions

- One process per `.nf` file in `modules/`
- Sub-workflows in `workflows/` using `take:` / `main:` / `emit:` DSL2 blocks
- Process aliases via `include { X as Y }` when reusing a process with different config
- Use `publishDir` overrides via `withName:` in config files with closures (`{ }`) for dynamic paths

### Adding a new module

1. Create a new `.nf` file in `modules/` (or an appropriate subdirectory)
2. Define a single process per file
3. Add resource configuration in the relevant config file (`assembly.config` or `analysis.config`)
4. Include the module in the appropriate pipeline or workflow

### Adding a new analysis block

1. Create a new workflow file in `workflows/`
2. Use the standard `take:` / `main:` / `emit:` structure
3. Add an enable flag in `analysis.config` (e.g., `my_block_enable = true`)
4. Wire it into `analysis.nf` guarded by the enable flag

## Testing

Before submitting a pull request, verify that your changes work by running the relevant pipeline on a small test dataset:

```bash
# Assembly pipeline
nextflow run assembly.nf -c assembly.config -profile test

# Analysis pipeline
nextflow run analysis.nf -c analysis.config --internal results/assembly
```

## Pull Requests

- Keep pull requests focused on a single change
- Provide a clear description of what the change does and why
- Reference any related issues
- Ensure the pipeline runs without errors on test data

## Reporting Issues

When reporting issues, please include:

- Nextflow version (`nextflow -version`)
- Docker version (`docker --version`)
- The command you ran
- The error message or unexpected behavior
- Relevant log files (`.nextflow.log`)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

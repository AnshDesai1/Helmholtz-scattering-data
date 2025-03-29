# NN-Born

This repository contains a finite element codebase for generating synthetic wave scattering data and solving both the forward and inverse Helmholtz problems in two dimensions. It was used in support of the experiments in the preprint:

**"A Neural Network Enhanced Born Approximation for Inverse Scattering"**
https://arxiv.org/abs/2503.01596

---

## Overview

The main script, `circles.py`, simulates wave scattering from randomly generated circular inhomogeneities embedded in a homogeneous background. The forward problem is solved using the finite element method (FEM) via NGSolve, and synthetic far-field patterns are generated for a variety of incident fields. These far-field data are stored along with the ground-truth scatterer and its Born approximation reconstruction.

The code is modular and extensible to support various scatterer geometries, including ellipses, blobs, L-shapes, U-shapes, Shepp-Logan phantoms, and more--all of which are readily available in the `matrix_form.py` and `custom_mesh.py` files.

---

## Features

- Synthetic data generation for 2D inverse scattering problems
- Finite element discretization using NGSolve
- Parametric mesh generation for diverse scatterer geometries
- Solving the Helmholtz equation with complex refractive index contrasts
- Discretized Born operator for linearized inversion
- Parallel sample generation using Python's `concurrent.futures`
- Output written in HDF5 format for downstream ML tasks
- Optional plotting for visual verification of each sample

---

## File Structure

- `circles.py`: Main driver script for sample generation and dataset writing
- `fem.py`: FEM solver routines for the Helmholtz equation and far-field pattern computation. Born approximation discretization via quadrature.
- `matrix_form.py`: Generates synthetic refractive index images for various scatterer shapes (circle, ellipse, blob, etc.)
- `custom_mesh.py`: Builds NGSolve meshes for the same shapes using constructive geometry routines

---

## Output Format

All data is saved in an HDF5 file containing:

- `image`: Ground-truth contrast function (n(x) - 1) stored as flattened vectors
- `approx`: Reconstructed contrast via the Born approximation
- `farfield.real`, `farfield.imag`: Real and imaginary parts of far-field data (flattened over incident and observation directions)

---

# Basic usage:
python3 circles.py \
  --base-index=0 \             # Starting index in the dataset
  --n-sample=10 \              # Number of samples to generate in this batch
  --n-sample-max=20000 \       # Maximum total number of samples in the dataset
  --data-file=results.hdf5 \   # Output HDF5 file path
  --generate-plots \           # Optional: enable plot generation
  --plot-dir=plots/ \          # Optional: directory to save plots
  --debug                      # Optional: enable debug-level logging


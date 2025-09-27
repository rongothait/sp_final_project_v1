# Symetric Non-negative Matrix Factorization (SymNMF)

## Overview
This project implements the **Symetric Non-negative Matrix Factorizaion (SymNMF)** clustering algorithm using both **C** and **Python**.
It provides a **hybrid implementation**:
- the mathematical core is written in **C** for performance
- A **Python-C API module** (`symnmf_mod`) exposes the C functions to Python.
- Python Scripts handle I/O, orchestration and analysis.
- Includes **K-Means** algorithm in Python for comparison

## Installation
1. **Build C Extension Module** - using `setup.py` run: `python3 setup.py build_ext --inplace`
2. **Compile the C program** - run `make` (alternitavely run `gcc symnmf.c -o symnmf -lm`)
# Next Steps for ERHAM Python Translation

The foundational structure (`Pyerham.py`) has been established with the `ErhamState` class to handle variables from Fortran `COMMON` blocks and the CLI arguments configured using `argparse`.

To complete the translation of the ~3500-line `erhamv16g-r3.for` code, the following milestones should be addressed:

## Phase 1: Input/Output Parsing and Data Formatting
The first priority is getting the data safely into the state objects and formatting it back out.

- [ ] Implement **`read_input`** (`SUBROUTINE INPUT`)
  - Parse `FILEIN` formatted parameters.
  - Populate the `ErhamState` arrays (e.g., transitions, frequencies, weights).
  - Handle the assignment of molecular constraints (`MQ`, `N1`, `N2`, `NC`).
- [ ] Implement **`predict`** and **`order`** (`SUBROUTINE PREDIC`, `SUBROUTINE ORDER`)
  - Translate the logic that predicts spectra.
  - Replicate the format printing (including the JPL catalog file format logic added in Rev 3).

## Phase 2: Core Physics and Hamiltonian Matrix Generation
This phase involves the physical modeling part of the program (Effective rotational Hamiltonian).

- [ ] Implement **`HAMILT`** (Hamiltonian generation)
- [ ] Implement matrix utility subroutines:
  - **`INDMAT`**
  - **`DMAT`**
  - **`BMAT`**
  - **`ASYMRO`** (Asymmetric rotor logic)
  - **`EVENODD`**

## Phase 3: Iterative Fit and Derivatives
This phase implements the least-squares fitting against experimental data.

- [ ] Implement **`iterate_fit`** (`SUBROUTINE ITER`)
  - Manage the iteration loops for the fitting process.
- [ ] Implement Derivative logic:
  - **`DERIV`**
  - **`DERPAR`**
- [ ] Implement Statistical testing (`SUBROUTINE FTEST`)

## Phase 4: Linear Algebra & Least Squares Optimization
*Note: Fortran manually implemented SVD and Eigen solvers. In Python, these should ideally be replaced by `scipy.linalg` or `numpy.linalg` equivalents wherever possible to save time and reduce errors.*

- [ ] Replace or rewrite Least Squares Subroutines:
  - **`LEASQU1`** (Least squares solver)
  - **`LSVDB`**, **`LSVDF`**, **`LSVG2`**, **`DROTG`**, **`VHS12`** (SVD related logic - evaluate replacing with `scipy.linalg.svd`)
- [ ] Replace or rewrite Eigen Solvers:
  - **`SEIGCX`**, **`HOUSCX`**, **`REVCX`**, **`INVITR`**, **`BISECT`** (Matrix diagonalization and eigenvectors - evaluate replacing with `scipy.linalg.eigh` or similar routines).

## Notes and Guidelines
- Remember that Fortran is 1-indexed while Python arrays are 0-indexed. Take special care when translating array bounds and loop ranges.
- Pass the `ErhamState` instance sequentially to the translated functions to avoid global variable state problems.
- Continually cross-reference against `erhaminstrv16g-r3.txt` to ensure physical parameters are translated correctly.

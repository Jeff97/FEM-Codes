# FEM-Codes

_A chronological collection of Fortran finite-element experiments built around CMake, PETSc, and VTK, progressing from baseline examples to nonlinear three-dimensional block compression._

[![GitHub stars](https://img.shields.io/github/stars/Jeff97/FEM-Codes?style=social)](https://github.com/Jeff97/FEM-Codes/stargazers)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Upstream: PFEMFort](https://img.shields.io/badge/upstream-PFEMFort-6f42c1)](https://github.com/chennachaos/PFEMFort)

---

## 📋 Overview

The eight numbered directories preserve successive FEM development snapshots. Together they document a path from a baseline Fortran/PETSc codebase to more general iteration strategies, hyperelastic constitutive models, three-dimensional elements, and a block-compression example.

This repository is derived from Chennakesava Kadapa's MIT-licensed [PFEMFort project](https://github.com/chennachaos/PFEMFort)[^1]. The original authorship and license are retained; see [Attribution and license](#-attribution-and-license).

## 📈 Learning progression

```mermaid
flowchart LR
    accTitle: FEM Development Progression
    accDescr: Eight numbered code snapshots progress from the baseline FEM implementation through general iteration and hyperelastic models to three-dimensional block compression.

    baseline_fem["1 · Baseline FEM"] --> january_snapshot["2 · January snapshot"]
    january_snapshot --> general_approach["3 · General approach"]
    general_approach --> nonlinear_iteration["4 · Nonlinear iteration"]
    nonlinear_iteration --> neo_hookean["5 · Neo-Hookean model"]
    neo_hookean --> mooney_rivlin["6 · Mooney-Rivlin model"]
    mooney_rivlin --> three_dimensional["7 · Three-dimensional FEM"]
    three_dimensional --> block_compression["8 · Block compression"]

    classDef foundation fill:#f3f4f6,stroke:#6b7280,stroke-width:2px,color:#1f2937
    classDef nonlinear fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f
    classDef advanced fill:#dcfce7,stroke:#16a34a,stroke-width:2px,color:#14532d

    class baseline_fem,january_snapshot,general_approach foundation
    class nonlinear_iteration,neo_hookean,mooney_rivlin nonlinear
    class three_dimensional,block_compression advanced
```

## 📚 Repository map

| Snapshot | Main focus |
| --- | --- |
| [`1 FEM/`](1%20FEM/) | Baseline FEM and PETSc examples |
| [`2 FEM20220124/`](2%20FEM20220124/) | Early reorganized snapshot |
| [`3 FEM20220204GeneralApproach/`](3%20FEM20220204GeneralApproach/) | Generalized implementation structure |
| [`4 FEM20220225Iteration/`](4%20FEM20220225Iteration/) | Iterative nonlinear workflow |
| [`5 FEM20220306NeoHookean/`](5%20FEM20220306NeoHookean/) | Neo-Hookean material model |
| [`6 FEM20220320MooneyRivlin/`](6%20FEM20220320MooneyRivlin/) | Mooney–Rivlin material model |
| [`7 FEM20220323For3D/`](7%20FEM20220323For3D/) | Three-dimensional elements |
| [`8 FEM20220404BlockCompression/`](8%20FEM20220404BlockCompression/) | Hexahedral block-compression example |

Several snapshots also contain generated `build/` artifacts. Treat the numbered directories as historical development records rather than interchangeable releases.

## 🔧 Build a snapshot

### Toolchain represented in the project

| Requirement | Role |
| --- | --- |
| Fortran compiler | Builds the FEM sources |
| CMake 3.0 or newer | Configures the included projects |
| PETSc | Linear algebra and solver infrastructure |
| VTK | Result-file support in the included utilities |
| BLAS and LAPACK | Numerical linear algebra dependencies |

The committed `CMakeLists.txt` files contain Linux include, library, and install paths from the original development environment. Update those paths for your system before configuring a build.

```bash
git clone https://github.com/Jeff97/FEM-Codes.git
cd FEM-Codes
cmake -S "8 FEM20220404BlockCompression" -B "8 FEM20220404BlockCompression/build-local"
cmake --build "8 FEM20220404BlockCompression/build-local"
```

The final snapshot defines the `hexelasticityimplicit` target. Other snapshots expose different targets; inspect the selected directory's `CMakeLists.txt` before building.

## 🔍 Reuse notes

- Start with one numbered snapshot and keep its source, input, and build assumptions together
- Prefer a new local build directory instead of modifying the committed historical `build/` tree
- Verify PETSc, VTK, BLAS, and LAPACK paths before diagnosing compiler errors
- Review author headers and the upstream project before redistributing modified source files

## 🔐 Attribution and license

This repository retains the [MIT License](LICENSE) from [PFEMFort](https://github.com/chennachaos/PFEMFort), including the original copyright notice for `chennachaos`. Source-file author notices, including those naming Dr. Chennakesava Kadapa, must also be preserved.

If this collection is useful for your FEM studies, consider starring this repository and the upstream project.

[^1]: Chennakesava Kadapa. “PFEMFort: Parallel programming for FEM using FORTRAN and PETSc.” https://github.com/chennachaos/PFEMFort

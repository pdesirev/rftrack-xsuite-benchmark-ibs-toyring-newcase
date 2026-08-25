# Intrabeam Scattering Benchmark — RF-Track vs Xsuite

**Authors:** Paula Desire, Andrea Latina

**Acknowledgements:** F. Soubelet, S. Kostoglou, H. Bartosik, G. Iadarola

## Overview

This repository contains the **latest version of a benchmark between the Intrabeam Scattering (IBS) models implemented in RF-Track and Xsuite**.

The objective is to compare the predictions of the different IBS implementations under equivalent beam and lattice conditions, with particular emphasis on the evolution of the transverse and longitudinal emittances.

The IBS models considered are:

### Xsuite

* **Kinetic IBS model** developed by Zampetakis et al.
* **Analytical IBS model** based on the Nagaitsev formalism.

### RF-Track

* **Monte Carlo-based IBS model** developed by P. Desire et al.

## Improvements with Respect to the Previous Benchmark

This repository is an updated version of the previous RF-Track/Xsuite IBS benchmark.

One of the main differences concerns the choice of beam and simulation parameters. In the previous benchmark, the IBS kicks could lead to stronger deviations from a Gaussian beam distribution.

In this updated case, the parameters have been adjusted so that the resulting distributions remain closer to **Gaussian distributions**. This provides a more suitable comparison between the different implementations, particularly because the analytical and kinetic approaches considered in the benchmark rely on assumptions associated with Gaussian beam distributions.

With these conditions, **better agreement between the RF-Track and Xsuite predictions is observed**.

The updated benchmark also includes:

* Tracking over a **larger number of turns**.
* Comparison of transverse and longitudinal emittance evolution.
* Analysis of the **6D emittance**.
* Application of **sigma cuts** to the particle distributions.
* Equivalent simulations **without IBS**, providing a reference for the IBS-induced emittance growth.

## Simulation Workflow

### 1. Generate the Initial Bunch

Use:

```bash
python bunch_converter.py
```

`bunch_converter.py` generates the initial bunch used for the benchmark and provides compatible initial conditions for the different tracking implementations.

Using equivalent initial distributions is essential for a meaningful comparison between RF-Track and Xsuite.

### 2. Xsuite — Kinetic IBS Model

Run:

```bash
python XSuite_tracking.py
```

This simulation uses the **kinetic IBS model developed by Zampetakis et al.**

### 3. Xsuite — Analytical IBS Model

Run:

```bash
python Xsuite_analytical.py
```

This simulation uses the **analytical Nagaitsev IBS model** implemented in Xsuite.

### 4. RF-Track — Monte Carlo IBS Model

Run:

```bash
octave rft.m
```

This simulation uses the **Monte Carlo-based IBS implementation in RF-Track**.

## Reference Simulations Without IBS

Equivalent simulations **without Intrabeam Scattering** are also included.

These provide a baseline against which the additional emittance growth produced by IBS can be evaluated.

The general comparison is therefore:

| Simulation framework | IBS model                                  |
| -------------------- | ------------------------------------------ |
| Xsuite               | Kinetic model — Zampetakis et al.          |
| Xsuite               | Analytical model — Nagaitsev               |
| RF-Track             | Monte Carlo-based model — P. Desire et al. |
| Xsuite / RF-Track    | No IBS — reference cases                   |

## Analysis

The different plotting scripts included in the repository can be used to analyze the evolution of the beam distribution and emittances.

The main quantities considered are:

* Horizontal transverse emittance
* Vertical transverse emittance
* Longitudinal emittance
* 6D emittance
* Emittance growth as a function of turn number

### Sigma Cuts

**Sigma cuts can be applied to the particle distributions** during the analysis.

This is particularly useful when comparing the different IBS approaches because the Monte Carlo method can generate a larger population of particles in the tails of the distribution.

Results can therefore be studied both with and without sigma cuts, providing additional information about the influence of outliers on the calculated emittance growth.

## Benchmark Strategy

The overall workflow can be summarized as:

```text
                         bunch_converter.py
                                │
                                ▼
                    Common initial distribution
                                │
              ┌─────────────────┼─────────────────┐
              │                 │                 │
              ▼                 ▼                 ▼
       Xsuite Kinetic    Xsuite Analytical     RF-Track
       Zampetakis et al.     Nagaitsev        Monte Carlo
              │                 │                 │
              └─────────────────┼─────────────────┘
                                │
                                ▼
                         Tracking results
                                │
                                ▼
                       Emittance analysis
                                │
                 ┌──────────────┼──────────────┐
                 ▼              ▼              ▼
             Transverse     Longitudinal    6D emittance
             emittances      emittance
                                │
                                ▼
                     RF-Track / Xsuite
                        IBS comparison
```

## Simulation Tools

### Xsuite

The Xsuite simulations use the IBS implementations available within the Xsuite ecosystem.

For additional information, refer to the [Xsuite documentation](https://xsuite.readthedocs.io/).

### RF-Track

RF-Track is used for the Monte Carlo-based IBS simulations.

For additional information about RF-Track, its physics models, available elements, and interfaces, refer to the [RF-Track Reference Manual](https://gitlab.cern.ch/rf-track/rf-track-reference-manual/-/raw/master/RF_Track_reference_manual.pdf?ref_type=heads).

## Authors and Acknowledgements

**Authors**

* Paula Desire
* Andrea Latina

**Acknowledgements**

The authors would like to thank:

* F. Soubelet
* S. Kostoglou
* H. Bartosik
* G. Iadarola

for their contributions, discussions, and support related to this benchmark.

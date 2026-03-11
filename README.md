
# DVPhiGen

A lightweight Monte Carlo event generator for exclusive $\phi$ electroproduction

\[
e p \rightarrow e' p' \phi \rightarrow e' p' K^+ K^-
\]

The generator produces events in **LUND format**, suitable for use with **GEMC / CLAS12 simulation pipelines**.

The cross section model is based on the empirical parametrization described in the CLAS12 proposal:

https://www.jlab.org/exp_prog/proposals/12/PR12-12-007.pdf

Specifically, the differential cross section parameterization implemented here follows:

- Eq. (37) – Eq. (66)

These equations describe the kinematic dependence of the exclusive $\phi$ electroproduction cross section as a function of:

- \(Q^2\)
- \(W\)
- \(t\)
- \(\phi\)
- \(\cos{\theta_{H}}\)

---


# Installation

Compile using the provided script:

```
./compile.sh
```

This produces the executable:

```
./build/gen
```

---

# Running the Generator

Default run:

```
./build/gen
```

Default parameters:

```
  EB       Beam energy in GeV        (default: 10.6)
  Q2min    Minimum Q^2 in GeV^2      (default: 1.0)
  Q2max    Maximum Q^2 in GeV^2      (default: 8.0)
  Wmin     Minimum W in GeV          (default: 2.0)
  Wmax     Maximum W in GeV          (default: 4.0)
  NEvents  Number of accepted events (default: 10000000)
  outfile  Output LUND file name     (default: dvphi.lund)
```

Output file:

```
dvphi.lund
```

---

# Command Line Arguments

All arguments:

```
./build/gen [EB] [Q2min] [Q2max] [Wmin] [Wmax] [NEvents] [outfile]
```

| Argument | Meaning | Default |
|--------|--------|--------|
| EB | Beam energy (GeV) | 10.6 |
| Q2min | Minimum \(Q^2\) | 1.0 |
| Q2max | Maximum \(Q^2\) | 8.0 |
| Wmin | Minimum \(W\) | 2.0 |
| Wmax | Maximum \(W\) | 4.0 |
| NEvents | Number of accepted events | 100000 |
| outfile | Output file name | dvphi.lund |

Example:

```
./build/gen 10.6 2 9 2 4 100000 dvphi.lund
```

---

# Help

```
./build/gen --help
```

# Example Output

```
========== Generator configuration ==========
EB      = 10.6 GeV
Q2 range= [1, 8] GeV^2
W  range= [2, 4] GeV
NEvents = 100000
Output  = dvphi.lund
====================

Progress: 10.0%
Progress: 20.0%
...
Progress: 100.0%

Generated 100000 accepted events into dvphi.lund

```

---

# Reference

CLAS12 φ production proposal:

https://www.jlab.org/exp_prog/proposals/12/PR12-12-007.pdf
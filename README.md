# PeriDyn

[![Coverage](https://codecov.io/gh/ravinderbhattoo/PeriDyn/branch/master/graph/badge.svg)](https://codecov.io/gh/ravinderbhattoo/PeriDyn) [![Generic badge](https://img.shields.io/badge/docs-ghpages-blue.svg)](https://ravinderbhattoo.github.io/PeriDyn)


<img src="https://ravinderbhattoo.github.io/PeriDyn/assets/diskdamage.png" alt="logo" style="display: block;
  margin-left: auto;
  margin-right: auto;
  width: 50%;"/>


PeriDyn is a numerical simulation software designed to solve peridynamics problems. It is written in the Julia programming language, and offers a high level of flexibility and speed. PDBenchmark is built on top of the PeriDyn package, which provides a number of predefined material models and benchmark problems. This allows users to quickly set up and run simulations, and compare their results to established benchmarks.


## Instructions

### Installation
```
import Pkg
Pkg.add(url="https://github.com/ravinderbhattoo/PDMaterialPoints")
Pkg.add(url="https://github.com/ravinderbhattoo/PeriDyn")
```

### Documentation
[HTML](https://ravinderbhattoo.github.io/PeriDyn)

[PDF](https://ravinderbhattoo.github.io/files/PeriDyn.pdf)

## Material models
1. Bond based peridynamics model
2. Ordinary state based peridynamics model
3. User defined material model


## Contact models
1. Linear and nonlinear repulsive model
2. User defined contact model


# Videos

|  Description | Video  |
|---|---|
|Disk Impact|<a href="http://www.youtube.com/watch?feature=player_embedded&v=RUdVr0Yh1jc " target="_blank"><img src="http://img.youtube.com/vi/RUdVr0Yh1jc/0.jpg" alt="IMAGE ALT TEXT HERE" width="480" height="360" border="1" /></a>|
|Wave propagation|<a href="http://www.youtube.com/watch?feature=player_embedded&v=q1N0aAdFYEs " target="_blank"><img src="http://img.youtube.com/vi/q1N0aAdFYEs/0.jpg" alt="IMAGE ALT TEXT HERE" width="480" height="360" border="1" /></a>|
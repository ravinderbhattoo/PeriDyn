# PeriDyn
Peridynamics implementation in Julia.

## Instructions

### Installation
```
import Pkg
Pkg.add(url="https://github.com/ravinderbhattoo/PDMesh.git")
Pkg.add(url="https://github.com/ravinderbhattoo/PeriDyn.git")
```

### Documentation
[HTML (incomplete)](https://ravinderbhattoo.github.io/PeriDyn)
[PDF](https://ravinderbhattoo.github.io/files/PeriDyn.pdf)

## Material models
1. Bond based peridynamics model
2. Ordinary state based peridynamics model
3. User defined material model


## Contact models
1. Linear and nonlinear repulsive model
2. User defined contact model

# Examples

## Block Plate Collision
![plate_block](/resources/plate_block.gif)

## Fracture
![fracture](/resources/notch.gif)

## Collision
![collision](/resources/2blocks.gif)

# Videos

|  Description | Video  |
|---|---|
|Disk Impact|<iframe width="500" height="400" src="https://www.youtube.com/embed/RUdVr0Yh1jc" title="Peridynamics: Disk Impact (bottom)" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>|
|Wave propagation|<iframe width="500" height="400" src="https://www.youtube.com/embed/q1N0aAdFYEs" title="Peridynamics: Stress Wave Propagation due to Impact Loading" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>|
# Contact models
In PeriDyn, there is a comprehensive set of contact models for peridynamics simulation. These contact models define the physical interactions between material particles and enable the accurate representation of contact forces and repulsion effects. In this post, we will explore the various contact models available in PeriDyn and provide an overview of the associated functions.

The repulsion models offer different formulations for calculating the repulsive forces between particles based on their relative positions. They are designed to simulate contact forces and prevent particle overlap, providing flexibility in modeling different contact behaviors. Such models can be applied to prevent overlap of particles of same material block and also between particles of different material blocks. The repulsion models are implemented as subtypes of the [`RepulsionModel`](#PeriDyn.RepulsionModel) abstract type. The following repulsion models are defined in PeriDyn:

```@docs
RepulsionModel11
RepulsionModel12
```

Here, `RepulsionModel11` is a repulsion model that calculates the repulsive forces between particles of the same material block. `RepulsionModel12` is a repulsion model that calculates the repulsive forces between particles of different material blocks.

## Short Range Repulsion Models

The ShortRangeRepulsionModel implements bond-based peridynamics forces for short-range repulsion interactions. It incorporates repulsive forces to ensure particles do not overlap within a specified range.
```@docs
ShortRangeRepulsionModel
```

## Simple Spring Models

PeriDyn also provides simple spring-like contact models that can be used in peridynamics simulations. These models describe the contact behavior between particles using spring forces based on their relative positions. The available models include:

```@docs
NonLinearSpringRepulsionModel
LinearSpringRepulsionModel
```

## LJ Models

The LJRepulsionModel represents a contact model based on the repulsive part of Lennard-Jones potential.

```@docs
LJRepulsionModel
```

## Functions for Repulsion Models

PeriDyn provides several functions that are utilized for working with repulsion models. These functions facilitate the implementation and calculation of repulsion forces within the contact models.

```@docs
repulsion_force
short_range_repulsion!
collision_box
update_repulsive_neighs!
RepulsionModel11_gcal
RepulsionModel12_gcal
```


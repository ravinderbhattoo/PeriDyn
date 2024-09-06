
# Contact models

In PeriDyn, there is a comprehensive set of contact models for peridynamics simulation. These contact models define the physical interactions between material particles and enable the accurate representation of contact forces and repulsion effects. In this post, we will explore the various contact models available in PeriDyn and provide an overview of the associated functions.

The contact models offer different formulations for calculating the repulsive forces between particles based on their relative positions. They are designed to simulate contact forces and prevent particle overlap, providing flexibility in modeling different contact behaviors. Such models can be applied to prevent overlap of particles of same material block and also between particles of different material blocks. The contact models are implemented as subtypes of the [`ContactModel`](@ref) abstract type. The following contact models are defined in PeriDyn:

- [`ContactModel11`](@ref)
- [`ContactModel12`](@ref)

Here, [`ContactModel11`](@ref) is a contact model that calculates the repulsive forces between particles of the same material block. [`ContactModel12`](@ref) is a contact model that calculates the repulsive forces between particles of different material blocks.

## Short Range Contact Models

The [`ShortRangeContactModel`](@ref) implements bond-based peridynamics forces for short-range contact interactions. It incorporates repulsive forces to ensure particles do not overlap within a specified range.

## Simple Spring Models

PeriDyn also provides simple spring-like contact models that can be used in peridynamics simulations. These models describe the contact behavior between particles using spring forces based on their relative positions. The available models include:

- [`NonLinearSpringContactModel`](@ref)
- [`LinearSpringContactModel`](@ref)

## Functions for Contact Models

PeriDyn provides several functions that are utilized for working with contact models. These functions facilitate the implementation and calculation of repulsion forces within the contact models.

- [`get_contact_force_fn`](@ref)
- [`short_range_contact!`](@ref)
- [`collision_box`](@ref)
- [`update_contact_neighs!`](@ref)
- [`ContactModel11_gcal`](@ref)
- [`ContactModel12_gcal`](@ref)

## Key Considerations

- **Neighbour Search**: To efficiently compute contact forces, contact models need to determine which material points are in contact. This is typically done using a neighbour search algorithm. The package includes functions like [`update_contact_neighs!`](@ref) that handle neighbour list updates for contact interactions.
- **Collision Detection**: Before calculating contact forces, the code often performs collision detection to determine if material blocks are potentially overlapping. This can be observed in functions like [`collision_box`](@ref), which defines a bounding box around material points to check for overlap.







# Material Models

# Material models in PeriDyn

A peridynamics material model in PeriDyn is defined using two key components: a **general material model** ([`GeneralMaterial`](@ref)) and a **specific material model**.  This approach separates the fundamental geometric and discretization data of a material from its constitutive behaviour.

The **general material model** ([`GeneralMaterial`](@ref)) stores information that dictates how a material is structurally represented within the simulation environment. This encompasses:

- **Particle Data**: This includes the positions (`x`, `y`), velocities (`velocity`), and volumes (`volume`) of individual particles that constitute the material. The deformed positions (`y`) can be initialised independently, allowing the simulation to start from a pre-deformed state.
- **Material Type and Block ID**: Each particle is associated with a type (`type`) and a block ID (`bid`).  The `type` distinguishes different materials within a simulation containing multiple materials, while the `bid` identifies particles belonging to the same material block.
- **Discretization Parameters**: The general model stores the particle size (`particle_size`) and the horizon (`horizon`). The horizon determines the range of interaction for each particle, a key aspect of peridynamic theory.
- **Neighbour Information**: The `family` array stores the indices of neighbouring particles within a particle's horizon. The `intact` array, a boolean array, indicates whether a bond between two particles is intact or broken based on the deformation criteria.
- **Weighted Volume**: The `weighted_volume` (represented as 'm' in some source code) is crucial for calculating internal forces.  It represents the influence of a particle's volume within its neighbourhood.

The **specific material model** defines the constitutive relationship of the material—essentially, how it responds to deformation.  These models encapsulate the material's intrinsic behaviour, independent of its geometric representation. Some specific material models include: 

- **[`BondBasedSpecific`](@ref)**: This model dictates material behaviour solely based on bond stretch between particle pairs. It is suitable for materials where forces arise primarily from bond elongation or compression. Key parameters include bulk modulus, critical stretch (the point at which a bond breaks), and density. The [`BondBasedMaterial`](@ref) combines this specific model with the general information to create a usable material block.
- **[`OrdinaryStateBasedSpecific`](@ref)**: This model incorporates both the deformation of individual bonds and the collective deformation of a particle's neighbourhood. It accounts for the influence of surrounding particles on a particle's deformation, making it suitable for a broader range of material behaviours compared to the bond-based model. This model utilizes parameters like bulk modulus, shear modulus, critical stretch, and density. Similar to the bond-based model, [`OrdinaryStateBasedMaterial`](@ref) combines the specific model ([`OrdinaryStateBasedSpecific`](@ref)) with the general material information.
- **[`ElastoPlasticSolidSpecific`](@ref)**: This model captures the behaviour of materials that exhibit both elastic and plastic deformation. It considers parameters such as bulk modulus, shear modulus, critical stretch, density, yield stress, and a plastic failure criterion.  The [`ElastoPlasticSolidMaterial`](@ref) combines this specific model with the general material information to model a block of material.

Each specific material model defines a [`force_density_t_ij`](@ref) function. This function calculates the force density, representing the internal forces acting on a particle due to its interaction with its neighbours, based on the specific material model used. 

This separation of general and specific material models in PeriDyn allows for a modular and flexible approach to defining materials. Users can easily switch between different constitutive models while keeping the geometric and discretization aspects consistent.  This facilitates the exploration of various material behaviours within the same simulation setup, simplifying comparative analyses and model development.



# Details of material models

PeriDyn offers a wide range of material models specifically designed for peridynamics simulations. These material models play a crucial role in defining the physical properties of material particles and ensuring accurate representation of material behavior. In this documentation, we will delve into the various material models available in PeriDyn and provide an overview of their associated functions and parameters.

The material models in PeriDyn are implemented as subtypes of the [`PeridynamicsMaterial`](@ref) abstract type. Each material model belonging to [`PeridynamicsMaterial`](@ref) encompasses the following fields:

- `name`: The name assigned to the material model.
- `type`: The type of material elements associated with the model.
- `blockid`: The block ID of the material model.
- `general`: The general parameters of a peridynamics material model.
- `specific`: The specific parameters particular to the material model.

## General Parameters

The general parameters of a material model are defined using the [`GeneralMaterial`](@ref) type. The [`GeneralMaterial`](@ref) type encompasses the following parameters:

#### Arguments
- `y0::AbstractArray{Float64,2}`: Initial displaced position of the material points.
- `v0::AbstractArray{Float64,2}`: Initial velocity of the material points.
- `x::AbstractArray{Float64,2}`: Initial position of the material points.
- `volume::AbstractArray{Float64,1}`: Volume of the material points.
- `type::AbstractArray{Int64,1}`: Type of the material points.
- `horizon::Float64`: Horizon of the material.

#### Keyword Arguments
- `max_neigh::Int64` = 100: Maximum number of neighbors.
- `particle_size::Float64` = 0: Particle size of the material.
- `skip_bb::Bool` = false: Skip bond-based material.

## Specific Parameters

The specific parameters of a material model are defined based on the subtype of the [`SpecificMaterial`](@ref) type. PeriDyn offers several [`PeridynamicsMaterial`](@ref) models, each defined by its respective [`SpecificMaterial`](@ref). Let's explore some of these material models:

### SkipMaterial

The [`SkipMaterial`](@ref) model is employed to bypass the material model definition for a block. Consequently, no forces will be generated in the material for a given deformation. The [`SkipSpecific`](@ref) parameters required to define this material model include:

- `density`: Density of the material.

### BondBasedMaterial

The [`BondBasedMaterial`](@ref) model is used to define material properties in bond-based peridynamics simulations. It requires the following parameters to define the [`BondBasedSpecific`](@ref):

- `bulkmodulus`: Bulk modulus of the material.
- `criticalstretch`: Critical stretch of the material.
- `density`: Density of the material.
- `horizon`: Horizon of the material.

### OrdinaryStateBasedMaterial

The [`OrdinaryStateBasedMaterial`](@ref) model is utilized to define material properties in ordinary state-based peridynamics simulations. To define the [`OrdinaryStateBasedSpecific`](@ref), the following parameters are required:

- `bulkmodulus`: Bulk modulus of the material.
- `shearmodulus`: Shear modulus of the material.
- `criticalstretch`: Critical stretch of the material.
- `density`: Density of the material.
- `horizon`: Horizon of the material.

### ElastoPlasticSolidMaterial

The [`ElastoPlasticSolidMaterial`](@ref) model is employed to define material properties in elasto-plastic peridynamics simulations. The [`ElastoPlasticSolidSpecific`](@ref) is defined using the following parameters:

- `bulkmodulus`: Bulk modulus of the material.
- `shearmodulus`: Shear modulus of the material.
- `criticalstretch`: Critical stretch of the material.
- `density`: Density of the material.
- `horizon`: Horizon of the material.
- `yieldstress`: Yield stress of the material.
- `criteria`: Yield criteria of the material, such as [`VonMises`](@ref), [`DruckerPrager`](@ref), etc. The default is [`VonMises`](@ref).



# Example 

To create an ordinary state-based material block in PeriDyn, you need to define both the general and specific material properties and then combine them. Here's a breakdown based on the sources:

### Defining the General Material Model ([`GeneralMaterial`](@ref))

The general material model stores the structural and geometric information about your material block. Here's an example:

```julia
resolution = 1.0 
x1, v1, y1, vol1, type1 = unpack(create(Cuboid([-5 5; -5 5; -5 5]); 
resolution=resolution, type=1))
horizon1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, horizon1, max_neigh=200)
```

**Explanation:**

1.  **`resolution = 1.0`**: This line sets the discretization resolution for your material block. A lower resolution means larger particles and a coarser representation.
2.  **`x1, v1, y1, vol1, type1 = ...`**: This line defines the initial positions (`x1`), velocities (`v1`), deformed positions (`y1`), volumes (`vol1`), and type (`type1`) of the particles in the material block. The `create(Cuboid(...))` function is assumed to generate these arrays based on the specified cuboid dimensions and resolution.
3.  **`horizon1 = 3.0*resolution`**: This line sets the horizon for the material block.  In this example, it's three times the particle spacing (`resolution`).
4.  **`mat_gen1 = GeneralMaterial(...)`**: Finally, a `GeneralMaterial` object (`mat_gen1`) is created using the information defined above. This object stores the general material properties.

### Defining the Specific Material Model ([`OrdinaryStateBasedSpecific`](@ref))

The specific material model defines the constitutive behaviour using an `OrdinaryStateBasedSpecific` object. Here's an example:

```julia
K, G = 1.0, 1.0 # Bulk and Shear Modulus
density1 = 1.0
cstretch1 = 0.5 # Critical stretch
mat_spec1 = OrdinaryStateBasedSpecific([K], [cstretch1], [density1]) 
```

**Explanation:**

1.  **`K, G = 1.0, 1.0`**: This line defines the bulk modulus (`K`) and shear modulus (`G`) of the material. 
2.  **`density1 = 1.0`**: This sets the density of the material.
3.  **`cstretch1 = 0.5`**: This line specifies the critical stretch of the material, a key parameter that governs when bonds between particles break under deformation.
4.  **`mat_spec1 = OrdinaryStateBasedSpecific(...)`**: This line creates an [`OrdinaryStateBasedSpecific`](@ref) object, which encapsulates the specific material behaviour using the provided parameters.

### Creating the Material Block ([`OrdinaryStateBasedMaterial`](@ref))

Now that you have defined the general and specific material properties, you can create the [`OrdinaryStateBasedMaterial`](@ref) block:

```julia
block1 = OrdinaryStateBasedMaterial(mat_gen1, mat_spec1; name="block 1")
```

**Explanation:**

1.  **`block1 = OrdinaryStateBasedMaterial(...)`**: This line combines the general material model (`mat_gen1`) and the specific material model (`mat_spec1`) to create an [`OrdinaryStateBasedMaterial`](@ref) object named `block1`. This object represents a material block with the specified properties ready to be used in a PeriDyn simulation.

This combined approach in PeriDyn provides a structured way to define materials, separating the constitutive behaviour from the geometric representation, which allows for flexibility and modularity in defining different material models.

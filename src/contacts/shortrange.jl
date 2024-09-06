"""
This module contains the definition of the short range contact model.
"""

export ShortRangeContactModel


"""
    ShortRangeContactModel(K, horizon, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; mfactor=1.0, kwargs...)

Create a short range contact model with the given stifness `K` and horizon `horizon`.

# Arguments
- `K`: Stifness of the contact model.
- `horizon`: Horizon of the contact model.
- `mat1`: Material of the first body.
- `mat2`: Material of the second body.

# Keyword Arguments
- `mfactor`: Multiplication factor for the stifness. Default is `1.0`.
- `kwargs...`: Keyword arguments for the `LinearSpringContactModel` constructor.
"""
function ShortRangeContactModel(K, horizon, mat1::PeridynamicsMaterial,
    mat2::PeridynamicsMaterial; mfactor=1.0, kwargs...)
    stifness = 18K/(pi*horizon^4) * mfactor
    LinearSpringContactModel(stifness,
            mat1::PeridynamicsMaterial,
            mat2::PeridynamicsMaterial; kwargs...)
end


"""
    ShortRangeContactModel(K, horizon, mat1::PeridynamicsMaterial; mfactor=1.0, kwargs...)

Create a short range contact model with the given stifness `K` and horizon `horizon`.

# Arguments
- `K`: Stifness of the contact model.
- `horizon`: Horizon of the contact model.
- `mat1`: Material of the first body.

# Keyword Arguments
- `mfactor`: Multiplication factor for the stifness. Default is `1.0`.
- `kwargs...`: Keyword arguments for the `LinearSpringContactModel` constructor.
"""
function ShortRangeContactModel(K, horizon, mat1::PeridynamicsMaterial;
    mfactor=1.0, kwargs...)
    stifness = 18K/(pi*horizon^4) * mfactor
    LinearSpringContactModel(stifness, mat1::PeridynamicsMaterial; kwargs...)
end

using ArgCheck

abstract type WiringDiagramAlgebra{T} end

"""
    combine(algebra, cod, dom1, dom2, map1, map2, arg1, arg2)

Let T be the operad of typed wiring diagrams. The arguments

    - `dom1`
    - `dom2`
    - `dom3`
    - `map1`
    - `map2`

specify a morphism D: dom1, dom2 → dom3 in T.

         + --------------- +
         |       dom3      |
         |        ↓ id     |
    D := |       dom3      |
         | map1 ↗   ↖ map2 |
         |  dom1     dom2  |
         + --------------- +

The argument `algebra` specifies a T-algebra, and the arguments

    - `arg1`
    - `arg2`

are elements of the sets D1 and D2.

    D1 := algebra(dom1)
    D2 := algebra(dom2)

This function computes the following element.

   algebra(D)(arg1, arg2) ∈ algebra(dom3)
"""
combine(
    algebra::WiringDiagramAlgebra{T},
    dom1::Integer,
    dom2::Integer,
    dom3::Integer,
    map1::AbstractVector,
    map2::AbstractVector,
    arg1::T,
    arg2::T,
) where {T}

"""
    project(algebra, dom1, dom2, map2, arg1)

Let T be the operad of typed wiring diagrams. The arguments

    - `dom1`
    - `dom2`
    - `map2`

specify a morphism D: dom1 → dom2 in T.

         + ------- +
         | dom2    |
         |  ↓ map2 |
    D := | dom1    |
         |  ↑ id   |
         | dom1    |
         + ------- +

The argument `algebra` specifies a T-algebra, and the argument

    - `arg1`

is an element of the set D1.

    D1 := algebra(dom1).

This function computes the following element.

   algebra(D)(arg1) ∈ algebra(dom2)
"""
project(
    algebra::WiringDiagramAlgebra{T},
    dom1::Integer,
    dom2::Integer,
    map2::AbstractVector,
    arg1::T,
)

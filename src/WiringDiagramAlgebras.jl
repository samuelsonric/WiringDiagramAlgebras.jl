module WiringDiagramAlgebras

using AbstractTrees
using ArgCheck
using Base: oneto
using Catlab: UndirectedWiringDiagram, TypedRelationDiagram, nparts
using CliqueTrees
using CliqueTrees: EliminationAlgorithm, DEFAULT_ELIMINATION_ALGORITHM, cliquetree!, nov
using CliqueTrees.Utilities
using FillArrays
using Graphs
using SparseArrays

export ClusterTree, WiringDiagramAlgebra, MatrixAlgebra, apply

include("wiring_diagram_algebra.jl")
include("matrix_algebra.jl")
include("static_wiring_diagram.jl")
include("cluster_tree.jl")

end

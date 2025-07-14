struct ClusterTree{T}
    types::Vector{T}
    root::Int
    VN::CliqueTree{Int, Int}                               # node → vertex
    VB::BipartiteGraph{Int, Int, Vector{Int}, Vector{Int}} # box  → vertex
    BN::BipartiteGraph{Int, Int, Vector{Int}, Vector{Int}} # node → box
end

function ClusterTree(diagram::UndirectedWiringDiagram; kwargs...)
    V = nparts(diagram, :Junction)
    weights = Ones{Int}(V)
    return ClusterTree(weights, diagram; kwargs...)
end

function ClusterTree(weights::AbstractVector, diagram::UndirectedWiringDiagram; alg::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM)
    P1 = nparts(diagram, :Port)
    P2 = nparts(diagram, :OuterPort)
    B = nparts(diagram, :Box) + 1
    V = nparts(diagram, :Junction)
    
    static = StaticWiringDiagram(diagram)
    types = static.types; VB = static.graph; BV = reverse(VB)
    
    matrix = sparse(BV)
    VV = BipartiteGraph(matrix' * matrix)

    perm, VN = cliquetree(weights, VV; alg)
    invp = invperm(perm); N = length(VN)
    permute!(types, perm)

    for p in oneto(P1 + P2)
        j = VB.tgt[p]
        VB.tgt[p] = invp[j]
    end

    root = 0; v = minimum(neighbors(VB, B))
    
    for (j, bag) in enumerate(VN)
        if v in residual(bag)
            root = j
            break
        end
    end
    
    perm = cliquetree!(VN, root)
    invp = invperm(perm)
    permute!(types, perm)    

    for p in oneto(P1 + P2)
        j = VB.tgt[p]
        VB.tgt[p] = invp[j]
    end
    
    roots = Vector{Int}(undef, V)
    
    for (j, bag) in enumerate(VN), v in residual(bag)
        roots[v] = j
    end
    
    tgt = Vector{Int}(undef, B)
        
    for b in oneto(B)
        v = minimum(neighbors(VB, b))
        tgt[b] = roots[v]
    end

    root = tgt[B]
    NB = BipartiteGraph(N, B, B, oneto(B + 1), tgt); BN = reverse(NB)
    return ClusterTree(types, root, VN, VB, BN)
end

function apply(algebra::WiringDiagramAlgebra{A, T}, tree::ClusterTree{T}, args) where {A, T}
    types = tree.types
    root = tree.root
    VN = tree.VN
    VB = tree.VB
    BN = tree.BN
    
    W = treewidth(VN) + 1
    N = length(VN)
    B = nov(BN)
    V = nov(VB)
    
    num = 0; stack = Vector{Tuple{Int, A}}(undef, N)
    
    query = neighbors(VB, B)

    typ3 = Vector{T}(undef, W)
    map2 = Vector{Int}(undef, W)
    prj3 = Vector{Int}(undef, V)
    prj4 = Vector{Int}(undef, V)
    marker = Vector{Int}(undef, V)
    
    for v in oneto(V)
        marker[v] = 0
    end
    
    for n in oneto(N)
        num = apply_collect_impl!(algebra, types, VN, VB, BN, args,
            stack, query, marker, typ3, map2, prj3, prj4, root, num, n)
    end
    
    # The variable `arg1` is an element of the set W1.
    #
    #     W1 := algebra(w1)
    #
    w1 = 0; arg1 = unit(algebra)::A
    
    while ispositive(num)
        # The variable `arg2` is an element of the set W2.
        #
        #     W2 := algebra(w2)
        #
        n, arg2 = stack[num]; num -= 1
        w1, arg1 = apply_combine_impl!(algebra, types, query, typ3, arg1, arg2, root, w1, n)
    end

    return arg1
end

function apply_collect_impl!(
        algebra::WiringDiagramAlgebra{A, T},
        types::AbstractVector{T},
        VN::CliqueTree{Int, Int},
        VB::BipartiteGraph{Int, Int},
        BN::BipartiteGraph{Int, Int},
        args,
        stack::AbstractVector{Tuple{Int, A}},
        query::AbstractVector{Int},
        marker::AbstractVector{Int},
        typ3::AbstractVector{T},
        map2::AbstractVector{Int},
        prj3::AbstractVector{Int},
        prj4::AbstractVector{Int},
        root::Int,
        num::Int,
        n::Int,
    ) where {A, T}
    B = nv(VB)
    
    # The variables
    #
    #     - `w4`
    #     - `inj4`
    #     - `prj4`
    #
    # specify a commutative diagram D.
    #
    #          + ------------- +
    #          |     inj4      |
    #          |    w4 → V     |
    #     D := | id ↓  ↙ prj4  |
    #          |    w4         |
    #          + ------------- +
    #
    w4 = 0; inj4 = VN[n]

    for v in inj4
        w4 += 1; prj4[v] = w4
    end

    # The variable `arg1` is an element of the set W1.
    #
    #     W1 := algebra(w1)
    #
    w1 = w3 = 0; arg1 = unit(algebra)::A

    for b in neighbors(BN, n)
        if b < B
            inj2 = neighbors(VB, b); arg2 = args[b]
            
            w1, arg1 = apply_collect_combine_impl!(algebra, types,
                marker, typ3, map2, inj2, prj3, arg1, arg2, w1, n)
        end
    end

    while ispositive(num)
        nn, arg2 = stack[num]

        n != parentindex(VN, nn) && break

        num -= 1; inj2 = separator(VN, nn)
        
        w1, arg1 = apply_collect_combine_impl!(algebra, types,
            marker, typ3, map2, inj2, prj3, arg1, arg2, w1, n)
    end

    # The variables
    #
    #     - `w1`
    #     - `w2`
    #     - `map2`
    #     - `inj2`
    #     - `prj3`
    #
    # specify a commutative diagram D.
    #
    #          + -------------- +
    #          |        inj2    |
    #          |      w2 → V    |
    #     D := | map2 ↓  ↙ prj3 |
    #          |      w1        |
    #          + -------------- +
    #
    w2 = 0

    if n == root
        inj2 = query
    else
        inj2 = separator(VN, n)
    end

    for v in inj2
        @assert marker[v] == n
        w2 += 1; map2[w2] = prj3[v]
    end

    # Let T be the operad of typed wiring diagrams. The variables
    #
    #    - `w1`
    #    - `w2`
    #    - `map2`
    #
    # specify a morphism D: w1 → w2 in T.
    #
    #          + ------ +
    #          | w2     |
    #          | ↓ map2 |
    #     D := | w1     |
    #          | ↑ id   |
    #          | w1     |
    #           + ------ +
    #
    # The variables
    #
    #    - `arg1`
    #    - `arg2`
    #
    # are elements of the sets W1 and W2.
    #
    #     W1 := algebra(w1)
    #     W2 := algebra(w2)
    #
    # They are related by the following equation.
    #
    #     arg2 = algebra(D)(arg1)
    #
    arg2 = project(algebra, w1, w2, map2, arg1, typ3)::A
    num += 1; stack[num] = (n, arg2)
    return num
end

function apply_collect_combine_impl!(
        algebra::WiringDiagramAlgebra{A, T},
        types::AbstractVector{T},
        marker::AbstractVector{Int},
        typ3::AbstractVector{T},
        map2::AbstractVector{Int},
        inj2::AbstractVector{Int},
        prj3::AbstractVector{Int},
        arg1::A,
        arg2::A,
        w1::Int,
        n::Int,
    ) where {A, T}
    # The variables
    #
    #     - `w2`
    #     - `w3`
    #     - `map2`
    #     - `inj2`
    #     - `prj3`
    #
    # specify a commutative diagram D.
    # 
    #          + -------------- +
    #          |        prj3    |
    #          |      w3 ← V    |
    #     D := | map2 ↑  ↗ inj2 |
    #          |      w2        |
    #          + -------------- +
    #
    # The variable `arg2` is an element of the set W2.
    #
    #     W2 := algebra(w2)
    #
    w2 = 0; w3 = w1

    for v in inj2
        if marker[v] < n
            marker[v] = n
            w3 += 1; prj3[v] = w3; typ3[w3] = types[v]
        end

        w2 += 1; map2[w2] = prj3[v]
    end

    # Let T be the operad of typed wiring diagrams. The variables
    #
    #    - `w1`
    #    - `w2`
    #    - `w3`
    #    - `map1`
    #    - `map2`
    #
    # specify a morphim D: w1, w2 → w3 in T.
    #
    #          + --------------- +
    #          |        w3       |
    #          |        ↓ id     |
    #     D := |        w3       |
    #          | map1 ↗   ↖ map2 |
    #          |    w1      w2   |
    #          + --------------- +
    #
    # The variables
    #
    #    - `arg1`
    #    - `arg2`
    #    - `arg3`
    #
    # are elements of the sets W1, W2, and W3.
    #
    #     W1 := algebra(w1)
    #     W2 := algebra(w2)
    #     W3 := algebra(w3)
    #
    # They are related by the following equation.
    #
    #     arg3 = algebra(D)(arg1, arg2)
    #
    map1 = oneto(w1)
    arg3 = combine(algebra, w1, w2, w3, map1, map2, arg1, arg2, typ3)::A
    return w3, arg3
end

function apply_combine_impl!(
        algebra::WiringDiagramAlgebra{A, T},
        types::AbstractVector{T},
        query::AbstractVector{Int},
        typ3::AbstractVector{T},
        arg1::A,
        arg2::A,
        root::Int,
        w1::Int,
        n::Int,
    ) where {A, T}
    
    # Let T be the operad of typed wiring diagrams. The variables
    #
    #    - `w1`
    #    - `w2`
    #    - `w3`
    #    - `map1`
    #    - `map2`
    #
    # specify a morphism w1, w2 → w3 in T.
    #
    #          + --------------- +
    #          |        w3       |
    #          |        ↓ id     |
    #     D := |        w3       |
    #          | map1 ↗   ↖ map2 |
    #          |    w1      w2   |
    #          + --------------- +
    #
    # The variables
    #
    #    - `arg1`
    #    - `arg2`
    #    - `arg3`
    #
    # are elements of the sets W1, W2, and W3.
    #
    #     W1 := algebra(w1)
    #     W2 := algebra(w2)
    #     W3 := algebra(w3)
    #
    # They are related by the following equation.
    #
    #     arg3 = algebra(D)(arg1, arg2)
    #
    w2 = 0

    if n == root
        for v in query
            w2 += 1; typ3[w2] = types[v]
        end
    end
    
    map1 = oneto(w1)
    map2 = oneto(w2)
    w3 = max(w1, w2)
    
    arg3 = combine(algebra, w1, w2, w3, map1, map2, arg1, arg2, typ3)::A
    return w3, arg3
end

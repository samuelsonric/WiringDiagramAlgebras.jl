struct StaticWiringDiagram{T}
    types::Vector{T}
    graph::BipartiteGraph{Int, Int, Vector{Int}, Vector{Int}}
end

function StaticWiringDiagram(diagram::StaticWiringDiagram)
    types = copy(diagram.types)
    graph = copy(diagram.graph)
    return StaticWiringDiagram(types, graph)
end

function StaticWiringDiagram(diagram::UndirectedWiringDiagram)
    V = nparts(diagram, :Junction); types = fill(nothing, V)
    return StaticWiringDiagram(types, diagram)
end

function StaticWiringDiagram(diagram::TypedRelationDiagram{T}) where {T}
    V = nparts(diagram, :Junction); types = Vector{T}(undef, V)

    for v in oneto(V)
        types[v] = diagram[v, :junction_type]
    end

    return StaticWiringDiagram(types, diagram)
end

function StaticWiringDiagram(types::AbstractVector{T}, diagram::UndirectedWiringDiagram) where {T}
    P1 = nparts(diagram, :Port)
    P2 = nparts(diagram, :OuterPort)
    B = nparts(diagram, :Box) + 1
    V = nparts(diagram, :Junction)

    ptr = Vector{Int}(undef, B + 1)
    tgt = Vector{Int}(undef, P1 + P2)

    ptr[1] = ptr[2] = 1

    for b in oneto(B - 1)
        ptr[b + 2] = 0
    end

    for p in oneto(P1)
        b = diagram[p, :box]
        ptr[b + 2] += 1
    end

    for b in oneto(B - 1)
        ptr[b + 2] += ptr[b + 1]
    end

    for p in oneto(P1)
        j = diagram[p, :junction]
        b = diagram[p, :box]
        tgt[ptr[b + 1]] = j; ptr[b + 1] += 1
    end

    for p in oneto(P2)
        j = diagram[p, :outer_junction]
        tgt[ptr[B + 1]] = j; ptr[B + 1] += 1
    end

    graph = BipartiteGraph(V, B, P1 + P2, ptr, tgt)
    return StaticWiringDiagram(types, graph)
end

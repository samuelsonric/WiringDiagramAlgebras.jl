struct StaticWiringDiagram
    graph::BipartiteGraph{Int, Int, Vector{Int}, Vector{Int}}
end

function StaticWiringDiagram(diagram::StaticWiringDiagram)
    graph = copy(diagram.graph)
    return StaticWiringDiagram(graph)
end

function StaticWiringDiagram(diagram::UndirectedWiringDiagram)
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
    return StaticWiringDiagram(graph)
end

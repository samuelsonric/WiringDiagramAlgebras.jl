struct MatrixAlgebra{T} <: WiringDiagramAlgebra{Vector{T}, Int} end
    
function combine(
        algebra::MatrixAlgebra{T},
        dom1::Int,
        dom2::Int,
        dom3::Int,
        map1::AbstractVector{Int},
        map2::AbstractVector{Int},
        arg1::AbstractVector{T},
        arg2::AbstractVector{T},
        typ3::AbstractVector{Int},
    ) where {T}
    @argcheck dom1 <= length(map1)
    @argcheck dom2 <= length(map2)
    @argcheck dom2 <= length(typ3)
    
    inds1 = Vector{Int}(undef, dom1)
    inds2 = Vector{Int}(undef, dom2)
    inds3 = Vector{Int}(undef, dom3)
 
    for i3 in oneto(dom3)
        inds3[i3] = typ3[i3]
    end

    for i2 in oneto(dom2)
        i3 = map2[i2]
        inds2[i2] = typ3[i3]
    end

    for i1 in oneto(dom1)
        i3 = map1[i1]
        inds1[i1] = typ3[i3]
    end

    arr1 = reshape(arg1, inds1...)
    arr2 = reshape(arg2, inds2...)
    arr3 = Array{T}(undef, inds3...)
    
    for inds3 in CartesianIndices(arr3)
        for i1 in oneto(dom1)
            i3 = map1[i1]
            inds1[i1] = inds3[i3]
        end
        
        for i2 in oneto(dom2)
            i3 = map2[i2]
            inds2[i2] = inds3[i3]
        end
        
        arr3[inds3] = arr1[inds1...] * arr2[inds2...]
    end
    
    arg3 = reshape(arr3, length(arr3))
    return arg3
end

function project(
        algebra::MatrixAlgebra{T},
        dom1::Int,
        dom2::Int,
        map2::AbstractVector{Int},
        arg1::AbstractVector{T},
        typ1::AbstractVector{Int},
    ) where {T}
    @argcheck dom1 <= length(typ1)
    @argcheck dom2 <= length(map2)

    inds1 = Vector{Int}(undef, dom1)
    inds2 = Vector{Int}(undef, dom2)
    
    for i1 in oneto(dom1)
        inds1[i1] = typ1[i1]
    end

    for i2 in oneto(dom2)
        i1 = map2[i2]
        inds2[i2] = typ1[i1]
    end
    
    arr1 = reshape(arg1, inds1...)
    arr2 = zeros(T, inds2...)
    
    for inds1 in CartesianIndices(arr1)
        for i2 in oneto(dom2)
            i1 = map2[i2]
            inds2[i2] = inds1[i1]
        end
        
        arr2[inds2...] += arr1[inds1]
    end
    
    arg2 = reshape(arr2, length(arr2))
    return arg2
end

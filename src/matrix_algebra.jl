struct MatrixAlgebra{T} <: WiringDiagramAlgebra{T} end
    
function combine(
        algebra::MatrixAlgebra{T},
        dom1::Integer,
        dom2::Integer,
        dom3::Integer,
        map1::AbstractVector,
        map2::AbstractVector,
        arg1::T,
        arg2::T,
    ) where {T}
    
    inds1 = Vector{Int}(undef, dom1)
    inds2 = Vector{Int}(undef, dom2)
    inds3 = Vector{Int}(undef, dom3)
    
    for i1 in oneto(dom1)
        i3 = map1[i1]
        inds3[i3] = size(arg1, i1)
    end
    
    for i2 in oneto(dom2)
        i3 = map2[i2]
        inds3[i3] = size(arg2, i2)
    end
    
    arg3 = T(undef, inds3...)
    
    for inds3 in CartesianIndices(arg3)
        for i1 in oneto(dom1)
            i3 = map1[i1]
            inds1[i1] = inds3[i3]
        end
        
        for i2 in oneto(dom2)
            i3 = map2[i2]
            inds2[i2] = inds3[i3]
        end
        
        arg3[inds3] = arg1[inds1...] * arg2[inds2...]
    end
    
    return arg3
end

function project(
        algebra::MatrixAlgebra{T},
        dom1::Integer,
        dom2::Integer,
        map2::AbstractVector,
        arg1::T,
    ) where {E, T <: AbstractArray{E}}
    
    inds2 = Vector{Int}(undef, dom2)
    
    for i2 in oneto(dom2)
        i1 = map2[i2]
        inds2[i2] = size(arg1, i1)
    end
    
    arg2 = T(undef, inds2...); fill!(arg2, zero(E))
    
    for inds1 in CartesianIndices(arg1)
        for i2 in oneto(dom2)
            i1 = map2[i2]
            inds2[i2] = inds1[i1]
        end
        
        arg2[inds2...] += arg1[inds1]
    end
    
    return arg2
end

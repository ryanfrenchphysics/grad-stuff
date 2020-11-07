

function fromfunction(func::Function, shape, type)
    arr = Array{type}(undef, shape)
    for I in CartesianIndices(arr)
        arr[I] = func(Tuple(I)...)
    end
    return arr
end

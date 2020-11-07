module JuliaNumpyGradient
using PyCall

export juliagradient

#np = pyimport("numpy")

function juliagradient(mat, dk)
    np = pyimport("numpy")
    return np.gradient(mat,dk)
end


end # module

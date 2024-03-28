m1 = rand(100,100)
m2 = ones(100,100)
m3 = rand(3,3)
v3 = vec(m3)
rows,cols = size(m3)
rm3 = reshape(v3, rows,cols)

rm3 == m3
function sort_and_convert1(M1,n)
    #convert matrix to a vector
    vec1 = vec(M1)
    rows,cols = size(M1)
    #save the order of the vector
    indices = sortperm(vec1, rev = true)
    ones_vector = ones(Int8, length(vec1))
    ones_vector[1:n] .= 0
    ones_vector = ones_vector[indices]
    m_out = reshape(ones_vector, rows, cols)
    return m_out
end

function sort_and_convert2(M1,n)
    #convert matrix to a vector
    vec1 = vec(M1)
    indices = vec(CartesianIndices(M1))
    tuple_vec = Tuple{Float16, CartesianIndex{2}}[]
    m_out = ones(Int8, size(M1))
    for i in 1:length(vec1)
        t = (vec1[i], indices[i])
        push!(tuple_vec,t)
    end
    sort!(tuple_vec, alg = PartialQuickSort(n), rev = true)
    for i in 1:n
        m_out[tuple_vec[i][2]] = 0
    end
    #save the order of the vector
    
    return m_out
end

# Example vector
vector = [3, 1, 4, 1, 5, 9]

# Sort the vector and get the indices of the sorted elements
sorted_indices = sortperm(vector)
sort!(vector)
# Reconstruct the original vector using the sorted indices
original_vector = vector[sorted_indices]

# Display the original vector and the reconstructed vector
println("Original Vector:")
println(vector)
println("Reconstructed Vector:")
println(original_vector)

# Example vector
vector = [3, 1, 4, 1, 5, 9]

# Sort the vector and get the indices of the sorted elements
sorted_indices = sortperm(vector)

# Reconstruct the original vector using the sorted indices
original_vector = vector[sorted_indices]

# Display the original vector and the reconstructed vector
println("Original Vector:")
println(vector)
println("Reconstructed Vector:")
println(original_vector)

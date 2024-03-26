using Plots

color_ramp = cgrad([:darkgoldenrod, :lightyellow, :green], [0, 5, 10])
big_matrix = rand(100, 100).*40
pasture_plot = heatmap(big_matrix, c = color_ramp, clim = (0, 10), colorbar = :none, axis = nothing, aspect_ratio=1.0)


function increment_heights(heights::Matrix{T}, increment_function::Function, day::Int) where T
    increment = increment_function(day)
    new_heights = heights.+increment
    return new_heights
end

function instantaneous_plant_growth(x)
    L = 50
    k = 0.1
    x₀ = 45
    return (k*L*exp(-k*(x-x₀)))/((1 + exp(-k*(x-x₀)))^2)
end

plant_height_mat = ones(Float16, 10,10)

plant_height_mat = increment_heights(plant_height_mat, instantaneous_plant_growth, 40)

function find_available_cells(pasture_area::Matrix{T}, threshold::Float16) where T
    above_height_threshold = Int.(pasture_area .>= threshold)
    return above_height_threshold
end

my_threshold = Float16(0.5)
my_matrix = rand(2,2)
find_available_cells(my_matrix, my_threshold)

function cattle_access_list(amp::Bool, pasture_size::Tuple{Int, Int})
    cattle_access_list = []
    subpasture_dims = (Int(pasture_size[1]/5), Int(pasture_size[2]/2))
    for i in 1:100
        if amp
            pasture = make_amp_subpastures(subpasture_dims)
        else
            pasture = ones(Int8, pasture_size[1], pasture_size[2])
        end
        push!(cattle_access_list, pasture)
    end
    return cattle_access_list
end

function make_amp_subpastures(subpasture_dims::Tuple{Int,Int})
    amp_subpastures = Matrix{Int8}[]
    random_start_position = rand(1:10)
    for i in 1:10
        # Create a 200 row x 500 col matrix of zeros
        mat = zeros(Int8, subpasture_dims[1], subpasture_dims[2])
        #add the matrix to the list of subpastures
        push!(amp_subpastures, mat)
    end
    #replace the matrix at the random start position with a matrix of ones
    amp_subpastures[random_start_position] = ones(Int8, subpasture_dims[1], subpasture_dims[2])
    #stitch the subpastures together to form the full pasture
    full_pasture = stitch_subpastures(amp_subpastures)
    #This should return one 1000 x 1000 matrix
    return full_pasture
end

function split_subpastures(full_matrix::Matrix{T}, subpasture_dims::Tuple{Int, Int}) where T
    # Horizontally split the full matrix into two halves of 500 columns each
    
    first_half, second_half = full_matrix[:, 1:subpasture_dims[2]], full_matrix[:, subpasture_dims[2]+1:end]

    # Vertically split each half into smaller matrices
    #The indexing works as follows: ((i-1)*100+1):(i*100) will give us the indices for the ith 100x500 matrix
    # for example, if i = 2, then ((2-1)*100+1):(2*100) = 101:200
    #the single : at the end of the index is used to select all columns
    first_half_matrices = [first_half[((i-1)*subpasture_dims[1]+1):(i*subpasture_dims[1]), :] for i in 1:5]
    second_half_matrices = [second_half[((i-1)*subpasture_dims[1]+1):(i*subpasture_dims[1]), :] for i in 1:5]

    # Combine the smaller matrices into a list
    subpasture_list = vcat(first_half_matrices, second_half_matrices)

    return subpasture_list
end

function rotate_cattle(pastures::Matrix{T}) where T
    subpastures = split_subpastures(pastures, (2,5))
    new_subpastures = circshift(subpastures, 1)
    new_pasture = stitch_subpastures(new_subpastures)
    return new_pasture
end

function stitch_subpastures(amp_subpastures::Vector{Matrix{T}}) where T
    # Split the list of matrices into two equal halves
    first_half = amp_subpastures[1:5]
    second_half = amp_subpastures[6:10]

    # Vertically concatenate the matrices in each half
    first_half_stitched = vcat(first_half...)
    second_half_stitched = vcat(second_half...)

    # Horizontally concatenate the two halves
    full_matrix = hcat(first_half_stitched, second_half_stitched)

    return full_matrix
end

my_pasture = make_amp_subpastures((2,5))
rotate_cattle(my_pasture)
subpastures = split_subpastures(my_pasture, (2,5))
    new_subpastures = circshift(subpastures, 1)
    new_pasture = stitch_subpastures(new_subpastures)
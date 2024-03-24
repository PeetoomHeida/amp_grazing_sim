#######
#=
This is a program to simulate a grazing landscape under AMP and conventional grazing. The intent is
to compare the two systems primarily in terms of their spatial autocorrelation of grazing effects.
Initially these grazing effects will likely be modelled as plant height. In this simulation we will
models plant growth as a linear process, and grazing as a removal of height at a constant rate per grazing event.
The plan is to also incorporate some variation in terms of the spatial distribution preferred vegetation.
=#

#Initial parameters of the simulation
growing_period = 100 #days
rotational_frequecy = 5 #days
#preference ranges from 0 to 1, with 0 being complete avoidance and 1 being complete preference
#the landscape is 100000 m x 100000 m (10 km x 10 km)
#the pixel size is 1 m
landscape_edge_length = 100000 #m
patch_size = 2 #m

#the average plant patch size is variable


#=
Steps for the simulation
    1. Create a landscape matrix with patches of preferred vegetation
    2. Create a matrix to store the plant height.
    3. Define a function to simulate the instantaneous growth of plants. Sigmoidal shape
    4. Define grazing intensity (eg. it takes off 1/2 of available plant height)
    5. Create a matrix to store the spatial distribution of pastures and subpastures
    6. Decide how many cells need to be eaten per time step for animal nutrition 
    (likely best approximated by the number of cells per subpasture / rotational frequency) = n_eaten
    7. Set a height threshold for plants to be 'available to be eaten'
    8. Create a matrix that defines the area that animals have access to.
    9. 

    At each time step:
    For all cells that animals have access to, a n_eaten cells get eaten. 
    Cells get eaten according to their relative preference value.
    This means that we have to rank all cells available to animals and then defoliate the top n cells 
    where n = n_eaten. 

    need to think about how to create the subpasture matrix and how to code the rotations after 5 days. 
    need to think about how to code conventional pastures in the landscape. 

=#

"""
This function is used to scale a matrix up by inserting blank rows/columns in between the existing rows
then filling the blanks in with the neighbouring values. 
# Arguments 
- `patch_matrix::Matrix{T}`: The original matrix to be scaled. For our usage it represents the preference values
of forage. 
- `patch_size::Int`: The size of the patches in the landscape
- `landscape_size::Int`: The length of one size of the square landscape
"""
function scale_matrix(patch_matrix::Matrix{T}, landscape_matrix::Matrix{T}) where T
    scale_factor = size(landscape_matrix,1)/size(patch_matrix,1)
    # Calculate new dimensions
    new_rows = size(patch_matrix, 1) * scale_factor
    new_cols = size(patch_matrix, 2) * scale_factor
    
    # Create a new matrix with scaled dimensions
    scaled_matrix = Matrix{T}(undef, new_rows, new_cols)
    
    # Fill in the values by replicating each element
    for i in 1:new_rows
        for j in 1:new_cols
            scaled_matrix[i, j] = patch_matrix[div(i-1, scale_factor) + 1, div(j-1, scale_factor) + 1]
        end
    end
    
    return scaled_matrix
end
"""
This function is used to create a landscape matrix with patches of preferred vegetation.
# Arguments
- `landscape_size::Int`: The length of one side of the square landscape
- `patch_size::Int`: The size of the patches in the landscape. The preference values
are assumed to be randomly distributed with a mean of 0.5, and bounded between 0 and 1.

"""
function patch_preferences(landscape_size::Int, patch_size::Int)
    # Create a matrix of the landscape size
    landscape_matrix = zeros(Float64, landscape_size, landscape_size)
    
    # Create a patch matrix of the patch size
    patch_matrix = clamp.(randn(landscape_size/patch_size, landscape_size/patch_size) .+ 0.5, 0, 1)
    
    # Scale the patch matrix to the landscape size
    scaled_matrix = scale_matrix(patch_matrix, landscape_matrix)
    
    return scaled_matrix
end


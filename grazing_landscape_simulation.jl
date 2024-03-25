import Pkg
Pkg.activate(".")

begin
    using Distributions
    using Plots
    using StatsPlots
end

#######
#=
This is a program to simulate a grazing landscape under AMP and conventional grazing. The intent is
to compare the two systems primarily in terms of their spatial autocorrelation of grazing effects.
Initially these grazing effects will likely be modelled as plant height. In this simulation we will
models plant growth as a linear process, and grazing as a removal of height at a constant rate per grazing event.
The plan is to also incorporate some variation in terms of the spatial distribution preferred vegetation.
=#

#Initial parameters of the simulation
grazing_period = 100 #days
rotational_frequecy = 5 #days
#preference ranges from 0 to 1, with 0 being complete avoidance and 1 being complete preference
#the landscape is 100000 m x 100000 m (10 km x 10 km)
#the pixel size is 1 m
subpasture_size = (500,200) #m x m
subpasture_area = 500*200 #m^2
n_eaten = 500*200/5 #cells eaten per day
landscape_edge_length = 100000 #m
patch_size = 2 #m
minimum_height = 3 #cm, the minimum height of plants to be eaten
day_of_cattle_release = 10 #day of the year when cattle are released
"""
This function models the instantaneous growth of plants. It is based on the derivative of the logistic growth function.
    Each value of x corresponds to a day, with the y-value being the plant growth that day.
    Peak growth occurs on day 45. 

    # Logistic Growth function
    `L/(1 + exp(-k*(x-x₀)))`
    
    # Parameters
    - `x::Int`: The day of the growing_period
    - `L::Int`: The maximum plant height
    - `k::Float64`: The growth rate
    - `x₀::Int`: The day of the year when growth is maximized
    # Returns
    - `y::Float64`: The change in plant height on day x
"""
function logistic_growth_curve(x)
    L = 50
    k = 0.1
    x₀ = 45
    return L/(1 + exp(-k*(x-x₀)))
end

function instantaneous_plant_growth(x)
    L = 50
    k = 0.1
    x₀ = 45
    return (k*L*exp(-k*(x-x₀)))/((1 + exp(-k*(x-x₀)))^2)
end

#= x=0:1:100
plot(plant_growth_curve.(x), label="Plant growth curve", xlabel="Days", ylabel="Plant height")
#the average plant patch size is variable
plot(logistic_growth_curve.(x), label="Logistic growth curve", xlabel="Days", ylabel="Plant height")
plot(instantaneous_plant_growth.(x), label="Logistic growth curve derivative", xlabel="Days", ylabel="Plant height")
 =#

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


    heights -> cattle_have_access -> find_available_cells -> output * preference_values -> find_grazed_cells -> 
    output * heights -> make_image -> increment_heights -> repeat
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
    new_rows = Int(size(patch_matrix, 1) * scale_factor)
    new_cols = Int(size(patch_matrix, 2) * scale_factor)
    
    # Create a new matrix with scaled dimensions
    scaled_matrix = zeros(Float16, new_rows, new_cols)
    
    # Fill in the values by replicating each element
    for i in 1:new_rows
        for j in 1:new_cols
            scaled_matrix[i, j] = patch_matrix[Int(div(i-1, scale_factor) + 1), Int(div(j-1, scale_factor) + 1)]
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
    landscape_matrix = zeros(Float16, landscape_size, landscape_size)
    
    # Create a patch matrix of the patch size
    patch_matrix = Float16.(clamp.(randn(Int(landscape_size/patch_size), Int(landscape_size/patch_size)) .+ 0.5, 0, 1))
    
    # Scale the patch matrix to the landscape size
    scaled_matrix = scale_matrix(patch_matrix, landscape_matrix)
    
    return scaled_matrix
end

"""
This function is used to create a matrix to identify which cells in a pasture are grazed.
    We assume that the cows will eat the n_eaten cells with the highest preference values.

"""
function find_grazed_cells(M::Array{Matrix{T}}, n_eaten::Int) where T
    ungrazed_cells = Matrix{Int8}[]
    for i in 1:length(M)
        #Negate the matrix (all values become negative), the find the nth smallest value (furthest from zero)
        #Then negate that value to get the threshold for grazing.
        threshold = -nthsmallest(-M[i], n_eaten)
        #Create a matrix of the same size as M, where each cell is 1 if the value in M is less than the threshold
        #This allows us to directly multiply the resulting matrix with the plant_height matrix to calculate the new height values
        escapes_grazing = (Int.(M[i] .< threshold))
        push!(grazed_cells, escapes_grazing)
    end
   
    return ungrazed_cells
end

function find_available_cells(M::Array{Matrix{T}}, threshold::Float64) where T
    available_cells = Matrix{Int8}[]
    for i in 1:length(M)
        M[i] = Int.(M[i] .> threshold)
        push!(available_cells, M[i])
    end
    return available_cells
end

function increment_heights(heights::Matrix{T}, increment_function::Function, day::Int) where T
    increment = increment_function(day)
    new_heights = heights.+increment
    return new_heights
end

function make_image(list_of_heights::Array{Matrix{T},1}, day::Int) where T
    nrows, ncols = 10, 10
    # Check if the number of matrices matches the grid size
if length(list_of_heights) != nrows * ncols
    error("Number of matrices (", length(list_of_heights), ") doesn't match grid size (", nrows * ncols, ")")
  end
  
  # Function to stack matrices horizontally
  function hstack_matrices(mats)
    return hcat(mats...)
  end
  
  # Function to stack matrices vertically
  function vstack_matrices(mats)
    big_mat = vcat(mats[1:ncols])
    for i in 2:nrows
      big_mat = vcat(big_mat, vcat(mats[(i*ncols + 1):((i+1)*ncols)]))
    end
    return big_mat
  end
  
  # Create a new empty matrix to store the stitched data
  big_matrix = zeros(size(first(list_of_heights))[1], size(first(list_of_heights))[2] * ncols)
  
  # Loop through rows and stack matrices horizontally
  for i in 1:nrows
    starting_index = (i-1)*ncols + 1
    ending_index = i*ncols
    row_matrices = list_of_heights[starting_index:ending_index]
    big_matrix[i, :] = hstack_matrices(row_matrices)
  end
  
  # Stack the rows of the big_matrix vertically
  big_matrix = vstack_matrices(big_matrix)
  color_ramp = cgrad([:brown, :yellow, :green], [0, 0.5, 1])
  # Create a plot with appropriate dimensions
  plot(big_matrix, aspect = ncols/nrows, c=color_ramp, clim=(0, 1))  # Adjust clim for your value range
  
  # Customize the plot (optional)
  xlabel("Columns")
  ylabel("Rows")
  title("Stitched Matrix (10x10)")
  
  # Display the plot
  jpeg("day_$(day).jpg")
end

cattle_access_continuous = ones(Int8, 1000,1000)

function populate_preference_list()
    preference_list = []
    for i in 1:100
        mat1 = patch_preferences(1000, 10)
        mat2 = patch_preferences(1000, 125)
        mat3 = mat1 .* mat2
        push!(preference_list, mat3)
    end
    return preference_list
end

"""
used to stitch the subpastures together to form the full pasture. Called in both the make_amp_subpastures function
and the rotate_cattle function.
"""
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
function plant_height_list()
    plant_heights = []
    for i in 1:100
        mat = zeros(Float16, 1000, 1000)
        push!(plant_heights, mat)
    end
    return plant_heights
end

function cattle_access_list(amp::Bool)
    cattle_access_list = []
    starting_subpastures = start_positions()
    for i in 1:100
        if amp
            mat = make_amp_subpastures(starting_subpastures[i])
        else
            mat = ones(Int8, 1000, 1000)
        end
        push!(cattle_access_list, mat)
    end
    return cattle_access_list
end
"""
This function is used to initialize the subpastures for the AMP system. It creates a list of 10 subpastures. 
    The pasture initially used for grazing is set by the random_start_position argument.
"""
function make_amp_subpastures(random_start_position::Int)
    amp_pastures = Matrix{Int8}[]
    for i in 1:10
        mat = zeros(Int8, 500, 200)
        push!(amp_pastures, mat)
    end
    amp_pastures[random_start_position] = ones(Int8, 500, 200)
    full_pasture = stitch_subpastures(amp_pastures)
    return full_pasture
end

"""
This function is used to split the full pasture into a list of subpastures.
    This helps with the rotation, as we then circshift the list of subpastures and reconstruct the full pasture.

"""
function split_subpastures(full_matrix::Matrix{T}) where T
    # Horizontally split the full matrix into two halves
    first_half, second_half = full_matrix[:, 1:500], full_matrix[:, 501:end]

    # Vertically split each half into smaller matrices
    first_half_matrices = [first_half[(i-1)*100+1:i*100, :] for i in 1:5]
    second_half_matrices = [second_half[(i-1)*100+1:i*100, :] for i in 1:5]

    # Combine the smaller matrices into a list
    matrices = vcat(first_half_matrices, second_half_matrices)

    return matrices
end
"""
This function is used to rotate the subpasture that the cattle are grazing on. It takes 
the list of 100 pastures as an argument. For each pasture it splits it into a list of matrices
that represent the subpastures, then it uses circshift to move the 'grazed' pasture up 
one index in the list, then it stitches the pastures back together and returns the full pasture.
    """
function rotate_cattle(pastures::Array{Matrix{T}}) where T
    rotated_pastures = []
    for i in 1:100
        subpastures = split_subpastures(pastures[i])
        new_subpastures = circshift(subpastures, 1)
        new_pasture = stitch_subpastures(new_subpastures)
        rotated_pastures[i] = new_pasture
    end
end
"""
Generates a list of random start positions for the AMP subpastures
"""
function start_positions()
    return rand(1:10, 100)
end

function simulate_grazing()
    preference_list = populate_preference_list()
    plant_heights = plant_height_list()
    cattle_accesses = cattle_access_list(true)
    for i in 1:grazing_period
        plant_heights = increment_heights.(plant_heights, instantaneous_plant_growth, i)
        cattle_access = cattle_accesses[i]
        preference = preference_list[i]
        available_cells = find_available_cells(cattle_access, 3)
        available_preferences = available_cells .* preference
        grazed_cells = find_grazed_cells.(available_preferences, n_eaten)
        #plant_heights[i] = plant_heights[i] .* available_cells
        plant_heights[i] = plant_heights[i] .* grazed_cells
        make_image(plant_heights, i)
        cattle_access_list = rotate_cattle(cattle_access_list)
    end
end

#= function create_gif()
    run(`convert -delay 10 -loop 0 day_*.jpg grazing_simulation.gif`)
end =#

simulate_grazing()
#create_gif()
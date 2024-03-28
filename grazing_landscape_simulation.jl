import Pkg
Pkg.activate(".")
Pkg.instantiate()

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
#n_eaten = Int(500*200/5) #cells eaten per day
landscape_edge_length = 10 #m
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
            #Int(div(i-1, scale_factor) + 1) is used to get the row index of the patch matrix.
            #The Int() in needed to convert the result of div to an integer, as indexes are integers.
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
    #Create a distribution of preference values
    preference_distribution = TruncatedNormal(0.5, 0.3, 0, 1)
    # Create a matrix of the landscape size
    landscape_matrix = zeros(Float16, landscape_size, landscape_size)
    
    # Create a patch matrix of the patch size
    patch_matrix = Float16.(rand(preference_distribution, Int(landscape_size/patch_size), Int(landscape_size/patch_size)))
    
    # Scale the patch matrix to the landscape size
    scaled_matrix = scale_matrix(patch_matrix, landscape_matrix)
    
    return scaled_matrix
end

function nth_smallest(A::AbstractArray{T,N}, n::Integer) where {T,N}
    threshold = sort(vec(A); alg=PartialQuickSort(n))[n]
    return threshold
end


  
"""
This function is used to create a matrix to identify which cells in a pasture are grazed.
    We assume that the cows will eat the n_eaten cells with the highest preference values.
    Grazed cells get a value of 0, while ungrazed cells get a value of 1.
"""
function find_grazed_cells(grazing_area::Matrix{T}, n_eaten::Int) where T
    #=
    #Negate the matrix (all values become negative), the find the nth smallest value (furthest from zero)
    #Then negate that value to get the threshold for grazing.
    threshold = -nth_smallest(-grazing_area, n_eaten)
    #Create a matrix of the same size as M, where each cell is 1 if the value in M is less than the threshold
    #This allows us to directly multiply the resulting matrix with the plant_height matrix to calculate the new height values
    #Need to also bin the value as larger than 0, or else the whole area gets grazed and the plant height goes to zero. 
    
    #TODO for some reason the plant heights don't recover in this edge case. 
    escapes_grazing = (Int.(grazing_area .< threshold && grazing_area .> 0))
    =#
        #convert matrix to a vector
    vec1 = vec(grazing_area)
    indices = vec(CartesianIndices(grazing_area))
    tuple_vec = Tuple{Float16, CartesianIndex{2}}[]
    escapes_grazing = ones(Int8, size(grazing_area))
    for i in 1:length(vec1)
        t = (vec1[i], indices[i])
        push!(tuple_vec,t)
    end
    sort!(tuple_vec, alg = PartialQuickSort(n_eaten), rev = true)
    for i in 1:n_eaten
        escapes_grazing[tuple_vec[i][2]] = 0
    end
    #save the order of the vector
        
    return escapes_grazing
end
"""
This function identifies the plant height cells that are above a certain height
    threshold. It returns a matrix of the same size as the pasture area, with 1s in the cells
    that are above the threshold and 0s in the cells that are below the threshold.
"""
function find_available_cells(pasture_area::Matrix{T}, threshold::Float16) where T
        above_height_threshold = Int.(pasture_area .>= threshold)
        return above_height_threshold
end
"""
This function is used to increment the plant heights in the pasture. It takes the current plant heights
    an increment function, and the day of the simulation as arguments. The increment function is used to
    calculate the amount of growth that occurs on that day.
"""
function increment_heights(heights::Matrix{T}, increment_function::Function, day::Int) where T
    increment = increment_function(day)
    new_heights = heights.+increment
    return new_heights
end

function stitch_matrices(matrices::Vector{<:AbstractMatrix})
    # Check if all matrices have the same dimensions
    if all(all(size(m) == size(matrices[1])) for m in matrices)
        n = size(matrices[1])[1]  # Get dimension size from first matrix
        result = zeros(n * 10, n * 10)  # Pre-allocate result matrix
  
        for i in 1:10
            for j in 1:10
                start_row = (i - 1) * n + 1
                end_row = i * n
                start_col = (j - 1) * n + 1
                end_col = j * n
                result[start_row:end_row, start_col:end_col] = matrices[(i - 1) * 10 + j]
            end
          end
        return result
    else
        error("All matrices must have the same dimensions")
    end
end

"""
Creates a heat maps of plant heights for each day of the simulation. The heat maps are saved as png files.
    The function takes a list of matrices of plant heights, the day of the simulation, and a boolean value
    that determines whether the heat map should be of the full landscape or just one subpasture.
"""
function make_image_heights(list_of_heights::AbstractArray{Matrix{T},1}, day::Int; full_landscape::Bool = true) where T
    if full_landscape
        big_matrix = stitch_matrices(list_of_heights)
    else
        big_matrix = list_of_heights[30]
    end
  
    color_ramp = cgrad([:black, :goldenrod4, :olivedrab4, :green], [0,0.01, 0.5, 1])
  # Create a plot with appropriate dimensions
  
    pasture_plot = heatmap(big_matrix, c = color_ramp, clim = (0, 40), aspect_ratio = :equal, colorbar_title = "Plant Height (cm)");  # Adjust clim for your value range
  
  # Customize the plot (optional)
    day_padded = lpad(string(day), 3, '0')
    # Display the plot
    png(pasture_plot, "/home/isaac/Coding/grazing_simulation/images/day_$(day_padded).png")
end

function make_image_pastures(list_of_heights::AbstractArray{Matrix{T},1}, day::Int) where T
    big_matrix = stitch_matrices(list_of_heights)
  
    color_ramp = cgrad([:white, :black], [0, 1])
  # Create a plot with appropriate dimensions
  
  pasture_plot = heatmap(big_matrix, c = color_ramp, clim = (0, 1), colorbar = :none, axis = nothing, aspect_ratio = 1.0);  # Adjust clim for your value range
  
  # Customize the plot (optional)
  day_padded = lpad(string(day), 3, '0')
    # Display the plot
    png(pasture_plot, "/home/isaac/Coding/grazing_simulation/images/day_$(day_padded).png")
end

"""
Creates a 100-element list of matrices of plant preferences. Each matrix is 1000 x 1000.
    The matrices are created by multiplying two matrices of the same size, one with a patch size of 10 and the other with a patch size of 125.
    This is used to create a smoother landscape that hopefully allows better visulaization of the grazing patterns.
    The resulting matrix is then appended to the preference_list.

    This should only be called one time at the start of the simulation.
"""
function populate_preference_list(pasture_size::Int, patchsize1::Int, patchsize2::Int)
    preference_list = Array{Float16}[]
    for i in 1:100
        mat1 = patch_preferences(pasture_size, patchsize1)
        mat2 = patch_preferences(pasture_size, patchsize2)
        mat3 = sqrt.(mat1 .* mat2)
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
    second_half = reverse(amp_subpastures[6:10])

    # Vertically concatenate the matrices in each half
    first_half_stitched = vcat(first_half...)
    second_half_stitched = vcat(second_half...)

    # Horizontally concatenate the two halves
    full_matrix = hcat(first_half_stitched, second_half_stitched)

    return full_matrix
end
"""
Generates a list of length 100 of matrices of plant heights. Each matrix is 1000 x 1000.
At the start of the growing season we assume the plant height to be 3, so that cattle grazing can start on 
day 1. The plant height will be incremented by the growth function each day as we call increment_heights() 
"""
function plant_height_list(pasture_size::Tuple{Int, Int})
    plant_heights = Array{Float16, 2}[]
    for i in 1:100
        mat = ones(Float16, pasture_size[1], pasture_size[2]) .* 3
        push!(plant_heights, mat)
    end
    return plant_heights
end

"""
Creates a list of 100 pastures for the simulation. If `amp` is true, the pastures are created using the AMP system, where the 
    1000 x 1000 pasture is divided into 10 subpastures, and one of the subpastures is selected to be grazed.
"""
function cattle_access_list(amp::Bool, pasture_size::Tuple{Int, Int})
    cattle_access_list = Array{Int8,2}[]
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
"""
This function is used to initialize the subpastures for the AMP system. It creates a list of 10 subpastures. 
    Each time the function is called it generates a random starting position for the cattle to start in. 
"""
function make_amp_subpastures(subpasture_dims::Tuple{Int,Int})
    amp_subpastures = Array{Int8,2}[]
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

"""
This function is used to split the full pasture into a list of subpastures.
    This helps with the rotation, as we then circshift the list of subpastures and reconstruct the full pasture.
    The returns a list of 10 subpastures. It contains a reverse() on the 2nd half of the list
    This is because the subpastures are stitched together in a specific order, and the order needs to be reversed
    Our subpasture matrix should look like:

    `
     [1 10;
      2 9;
      3 8;
      4 7;
      5 6]. 
   `

    In order to have the list work with the circshift function we 
        need to reverse the order of the 2nd half of the list such that the list is returned as [1,2,3,4,5,6,7,8,9,10]

"""
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
    subpasture_list = vcat(first_half_matrices, reverse(second_half_matrices))

    return subpasture_list
end
"""
This function is used to rotate the subpasture that the cattle are grazing on. It takes 
the list of 100 pastures as an argument. For each pasture it splits it into a list of matrices
that represent the subpastures, then it uses circshift to move the 'grazed' pasture up 
one index in the list, then it stitches the pastures back together and returns the full pasture.
    """
function rotate_cattle(pastures::Matrix{T}) where T
    subpasture_dims = (Int(size(pastures,1)/5), Int(size(pastures,2)/2))
    #The returns a list of 10 subpastures. It contains a reverse() on the 2nd half of the list
    #This is because the subpastures are stitched together in a specific order, and the order needs to be reversed
    #Our subpasture matrix should look like [1,10;,2,9;3,8;4,7;5,6]. In order to have the list work with the circshift
    #function we need to reverse the order of the 2nd half of the list such that the list is returned as [1,2,3,4,5,6,7,8,9,10]
    subpastures = split_subpastures(pastures, subpasture_dims)
    #first_half, second_half = subpastures[1:5], subpastures[6:10]
    #This creates a circular motion in the pasture rotations
    new_subpastures = circshift(subpastures, 1)
    new_pasture = stitch_subpastures(new_subpastures)
    return new_pasture
end
"""
Generates a list of random start positions for the AMP subpastures
"""
function start_positions()
    return rand(1:10, 100)
end

function create_animation(amp_on::Bool)
    amp_string = amp_on ? "amp" : "conventional"
    run(`ffmpeg -framerate 3 -pattern_type glob -i './images/*.png' -c:v libx264 -pix_fmt yuv420p ./images/vegetation_$(amp_string).mp4
    `)

    run(`ffmpeg -framerate 10 -pattern_type glob -i './images/day_*.png' -vf "scale=800:-1" ./images/vegetation_$(amp_string).gif -y`
    )
end

function simulate_grazing()
    #subpasture_rows = 200
    pasture_size = 100
    subpasture_dims = (Int(pasture_size/5), Int(pasture_size/2))
    amp_on = true
    n_eaten = Int(round(subpasture_dims[1]*subpasture_dims[2]/5))
    preference_list = populate_preference_list(pasture_size, 20, 5)
    plant_heights = plant_height_list((pasture_size, pasture_size))
    cattle_accesses = cattle_access_list(amp_on, (pasture_size, pasture_size))
    for i in 1:100#grazing_period
        ##increment the plant heights for each pasture in the list
        plant_heights = increment_heights.(plant_heights, instantaneous_plant_growth, i)
        
        #print("max plant height is $(maximum(plant_heights[1]))\n")
        #multiply the cattle access matrix by the plant heights to get the accessible plant heights
        accesible_plants = [a .* b for (a,b) in zip(cattle_accesses, plant_heights)] 
        
        #identify all cells with a plant height >= 3
        prey_plants = find_available_cells.(accesible_plants, Float16(3))
        
        #print("available cells for day $(i): $(available_cells[1])\n")
        #multiply the available cells by the preference list to get the preference values for each cell
        #unavailable cells will have their preference values set to 0
        #available cells values will be preference * 1
        prey_plants_preferences = [a .* b for (a,b) in zip(prey_plants, preference_list)]
        
        #print("available_preferences: $(available_preferences[1])\n")
        ### TODO: change the n_eaten back into 20000
        grazed_cells = find_grazed_cells.(prey_plants_preferences, n_eaten)
        
        #print("grazed cells: $(grazed_cells[1])\n")
        
        #multiply the plant heights by the grazed cells to get the new plant heights
        #grazing sets the plant height to zero. 
        plant_heights = [a .* b for (a,b) in zip(plant_heights, grazed_cells)]
        #print("plant heights are now: $(plant_heights[1])\n")
        #To make images of plant height
        make_image_heights(plant_heights, i, full_landscape = false)

        #to make images of the rotation uncomment the line below
        #make_image_pastures(cattle_accesses, i)
        
        #print("saved image for day $(i)")
        
        if i % rotational_frequecy == 0 && amp_on
            cattle_accesses = rotate_cattle.(cattle_accesses)
        end
        
    end
    create_animation(amp_on)
end

#= function create_gif()
    run(`ffmpeg -f image2 -pattern_type glob -i '*.png' output_amp.gif
`)
end =#

simulate_grazing()
#create_gif()

# TODO Why do some fields never grow high?
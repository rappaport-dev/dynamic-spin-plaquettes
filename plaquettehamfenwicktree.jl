using ArgParse
using CSV
using DataFrames
using Random # We need this for rand()
using DataStructures

# Converts 2D (i,j) coordinates to a 1D index
function to_1d(i, j, L)
    return (j - 1) * L + i
end

# Converts a 1D index back to 2D (i,j) coordinates
function to_2d(index, L)
    i = mod1(index, L) # mod1 is perfect for this
    j = ((index - 1) ÷ L) + 1
    return (i, j)
end

function setup_model(L, p)

    # randomly assign how the hamiltonian looks
    # on 2x2 might be:
    # 1   5 
    # 1   1 
    # where due to periodic boundary conditions, the 5 will loop around

    hamiltonian_placements = [(rand() < p ? 1 : 5) for _ = 1:L, _ = 1:L]

    # create spin up and spin downs with equal probability 
    spins = [(rand() < 0.5 ? 1 : -1) for _ = 1:L, _ = 1:L]

    return (hamiltonian_placements, spins)

end

function calculate_total_energy(L, hamiltonian_placements, spins)

    total_energy = 0.0

    for j = 1:L, i = 1:L

        # local energy term 
        h_i = 0

        h_i = calculate_h_i(i, j, L, hamiltonian_placements, spins)

        total_energy += (1 - h_i)
    end

    return total_energy
end

function find_neighbors(i, j, L)
    # calculate neighbor indices with periodic boundaries
    up_coords = (mod1(i - 1, L), j)
    down_coords = (mod1(i + 1, L), j)
    left_coords = (i, mod1(j - 1, L))
    right_coords = (i, mod1(j + 1, L))

    # We can return a list (or tuple) of these coordinate tuples
    return [up_coords, down_coords, left_coords, right_coords]
end

function calculate_h_i(i, j, L, hamiltonian_placements, spins)

    h_i = 0.0 # Initialize h_i

    if hamiltonian_placements[i, j] == 1
        h_i = spins[i, j]
    else

        # spin at the site 
        h_i = spins[i, j]

        # loop over the 4 neighbor positions 
        for (neighbori, neighborj) in find_neighbors(i, j, L)
            h_i = h_i * spins[neighbori, neighborj]
        end

    end

    return h_i
end

function calculate_delta_E(i, j, L, hamiltonian_placements, spins)
    # calculates the total energy change after one spin flip 

    delta_E = 0.0 #initialize delta E 

    # 2*h_i term in change in energy
    delta_E += 2 * calculate_h_i(i, j, L, hamiltonian_placements, spins)

    # check up, down, left, right neighbors and check whether to add their term 
    for (neighbori, neighborj) in find_neighbors(i, j, L)
        if hamiltonian_placements[neighbori, neighborj] == 5
            delta_E +=
                2 * calculate_h_i(neighbori, neighborj, L, hamiltonian_placements, spins)
        end
    end

    return delta_E

end

function initialize_rates(L, hamiltonian_placements, spins, beta)
    total_rate = 0.0

    # matrix containing flip rates at each site 
    # flip rate is now a fenwick tree
    flip_rates = FenwickTree{Float64}(L * L)

    for j = 1:L, i = 1:L

        delta_E = calculate_delta_E(i, j, L, hamiltonian_placements, spins)

        w_i = 1 / (1+exp(beta * delta_E))

        total_rate += w_i
        inc!(flip_rates, to_1d(i, j, L), w_i)
    end

    return (flip_rates, total_rate)

end

function mcmc_step!(L, hamiltonian_placements, spins, beta, flip_rates, total_rate)
    """New version of mcmc step, it takes in a flip_rates and total_rates. """
    """
    Performs one rejection-free MCMC step using pre-calculated rates.
    Modifies spins and flip_rates in place.
    Returns the coordinates of the flipped spin (i, j), the time step delta_T,
    and the *updated* total_rate
    Now uses a FenwickTree as the datastructure 
    """

    # Handle frozen state
    if total_rate == 0.0
        # Return dummy coordinates, infinite time, and unchanged rate
        return ((nothing, nothing), Inf)
    end

    #calculate delta T, 
    delta_T = - log(rand()) / total_rate

    # calculate the chosen flipped spin with the algorithm mentioned in the supp
    spin_choice = rand()*total_rate

    # datastructures replaces old code with Fenwick Tree- finds which spin to flip in O(log(N)) time 
    index_1d = DataStructures.index_at_prefixsum(flip_rates, spin_choice)
    (flipped_i, flipped_j) = to_2d(index_1d, L)

    old_rate_k = flip_rates[index_1d] - (index_1d == 1 ? 0.0 : flip_rates[index_1d-1])

    # Check if a spin was actually found (should always happen if total_rate > 0)
    if flipped_i == 0
        # This case should ideally not be reached if total_rate > 0
        println("Warning: Spin selection in mcmc_step failed unexpectedly.")
        return ((nothing, nothing), delta_T, total_rate)
    end

    spins[flipped_i, flipped_j] *= -1
    # 1. Calculate the new rate for the flipped spin
    new_delta_E_k =
        calculate_delta_E(flipped_i, flipped_j, L, hamiltonian_placements, spins)
    new_rate_k = 1 / (1 + exp(beta * new_delta_E_k))

    # 2. Update the FenwickTree. This is an O(log N) operation.
    DataStructures.update!(flip_rates, index_k_1d, new_rate_k)

    # 3. Update the total_rate incrementally
    total_rate = total_rate - old_rate_k + new_rate_k

    # have to update the neighbor's spin flips too
    for (ni, nj) in find_neighbors(flipped_i, flipped_j, L)

        # 1. Get the 1D index for this neighbor
        index_j_1d = to_1d(ni, nj, L)

        # 2. Get the neighbor's OLD rate (using the O(log N) prefix sum trick)
        old_rate_j =
            flip_rates[index_j_1d] - (index_j_1d == 1 ? 0.0 : flip_rates[index_j_1d-1])

        # 3. Calculate the neighbor's NEW rate (this is unchanged)
        new_delta_E_j = calculate_delta_E(ni, nj, L, hamiltonian_placements, spins)
        new_rate_j = 1 / (1 + exp(beta * new_delta_E_j))

        # 4. Update the FenwickTree with the new rate
        DataStructures.update!(flip_rates, index_j_1d, new_rate_j)

        # 5. Update the total_rate incrementally
        total_rate = total_rate - old_rate_j + new_rate_j

    end

    return ((flipped_i, flipped_j), delta_T, total_rate)

end

function run_simulation(L, p, beta, T_max, snapshot_interval)
    """
    Runs a full simulation using the OPTIMIZED O(1) MCMC step.
    """
    hamiltonian_placements, spins = setup_model(L, p)

    current_time = 0.0
    next_snapshot_time = snapshot_interval
    snapshot_list = []

    # not sure why you would do this, but if T_max is 0 or less, return empty list 
    if T_max <= 0
        return snapshot_list
    end

    flip_rates, total_rate = initialize_rates(L, hamiltonian_placements, spins, beta)


    while current_time < T_max

        # perform a step, and get the spin coordinates and the change in time 
        ((flipped_i, flipped_j), delta_T, total_rate) =
            mcmc_step!(L, hamiltonian_placements, spins, beta, flip_rates, total_rate)

        # in this case, loop in mcmc_step failed 
        if flipped_i == 0
            print("Mcmc step did not flip a spin, simulation stopping")
            break
        end
        # i think this isn't possible anymore- would get error in FenwickTree?

        current_time += delta_T

        while current_time >= next_snapshot_time

            push!(snapshot_list, copy(spins))
            next_snapshot_time += snapshot_interval

        end

    end

    return snapshot_list
end




function calculate_correlation(snapshot_list, snapshot_interval, L, t_w)

    N_spins = L*L

    # need to cast index as int 
    t_w_index = Int(round(t_w / snapshot_interval))
    # also, make sure index is btwn 1 and the length of the lsit 
    t_w_index = clamp(t_w_index, 1, length(snapshot_list))


    spins_at_tw = snapshot_list[t_w_index]

    correlation_values = []

    # calculates all the intermediate values of the function C(t,tw)
    for t_prime_index = t_w_index:length(snapshot_list)

        spins_at_t_prime = snapshot_list[t_prime_index]

        # correlation function implementation 
        # .* is the element-wise product of two arrays. so just sum over the element-wise product of two arrays 
        correlation = sum(spins_at_tw .* spins_at_t_prime)

        # normalize by N_spins
        correlation /= N_spins

        push!(correlation_values, correlation)

    end

    return correlation_values

end

"""
Averages a list of curves (which may have different lengths).
It truncates all curves to the minimum length found and
then computes the element-wise average.
"""
function average_curves(list_of_curves)
    # Handle the case of an empty list
    if isempty(list_of_curves)
        return []
    end

    # 1. Find the number of runs
    num_runs = length(list_of_curves)

    # 2. Find the minimum length
    #    (This comprehension finds all lengths, then minimum finds the smallest)
    min_len = minimum(length(curve) for curve in list_of_curves)

    # 3. Initialize the accumulator array
    avg_curve = zeros(min_len)

    # 4. Loop through each curve from a single run
    for curve in list_of_curves
        # 5. Loop from 1 to the shortest length
        for i = 1:min_len
            # 6. Add this curve's contribution
            avg_curve[i] += curve[i]
        end
    end

    # 7. Normalize by the number of runs
    avg_curve /= num_runs

    # 8. Return the final averaged curve
    return avg_curve
end

function run_experiment(L, p, beta, T_max, snapshot_interval, t_w1, t_w2, num_runs)

    all_runs_c_tw1 = []
    all_runs_c_tw2 = []

    for i = 1:num_runs

        println("Running simulation $i of $num_runs")

        snapshot_list = run_simulation(L, p, beta, T_max, snapshot_interval)
        correlations_t_w1 = calculate_correlation(snapshot_list, snapshot_interval, L, t_w1)
        correlations_t_w2 = calculate_correlation(snapshot_list, snapshot_interval, L, t_w2)

        push!(all_runs_c_tw1, correlations_t_w1)
        push!(all_runs_c_tw2, correlations_t_w2)

    end

    #return averaged results 

    return (average_curves(all_runs_c_tw1), average_curves(all_runs_c_tw2))

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L", "-L"
        help = "Linear lattice size"
        arg_type = Int
        default = 16
        "--p", "-p"
        help = "Probability of 1-body term"
        arg_type = Float64
        default = 0.1
        "--beta", "-b"
        help = "Inverse temperature"
        arg_type = Float64
        default = 2.0
        "--T_max"
        help = "Total simulation time"
        arg_type = Float64
        default = 1000.0
        "--snapshot_interval"
        help = "Time interval between snapshots"
        arg_type = Float64
        default = 1.0
        "--t_w1"
        help = "First wait time for correlation"
        arg_type = Float64
        default = 10.0
        "--t_w2"
        help = "Second wait time for correlation"
        arg_type = Float64
        default = 100.0
        "--num_runs", "-n"
        help = "Number of runs for averaging"
        arg_type = Int
        default = 100
        "--output_file", "-o"
        help = "Path to save the output CSV file"
        arg_type = String
        default = "results.csv"
    end

    return parse_args(s)
end

function main()
    # 1. Get all parameters from the command line
    println("Parsing arguments...")
    parsed_args = parse_commandline()

    # 2. Run the full, averaged experiment
    println("Starting experiment...")
    (avg_c_tw1, avg_c_tw2) = run_experiment(
        parsed_args["L"],
        parsed_args["p"],
        parsed_args["beta"],
        parsed_args["T_max"],
        parsed_args["snapshot_interval"],
        parsed_args["t_w1"],
        parsed_args["t_w2"],
        parsed_args["num_runs"],
    )
    println("Experiment complete. Saving results...")

    # 3. Save the results to a CSV file

    # We need to create the time axis for the plot
    # C(t_w, t') is a function of (t' - t_w)
    # Both curves start at (t' - t_w) = 0

    # Find the shortest curve, as average_curves did
    min_len_1 = length(avg_c_tw1)
    min_len_2 = length(avg_c_tw2)

    # Create time axes. 
    # The time step is snapshot_interval, but the *index* is just 0, 1, 2...
    # Let's just save the data vs. the time difference.

    time_axis_1 = (0:(min_len_1-1)) .* parsed_args["snapshot_interval"]
    time_axis_2 = (0:(min_len_2-1)) .* parsed_args["snapshot_interval"]

    # Create a DataFrame (like a spreadsheet)
    # We pad the shorter list with 'missing' so they can be in one file
    max_len = max(min_len_1, min_len_2)

    df = DataFrame(
        delta_t_1 = vcat(time_axis_1, fill(missing, max_len - min_len_1)),
        C_tw1 = vcat(avg_c_tw1, fill(missing, max_len - min_len_1)),
        delta_t_2 = vcat(time_axis_2, fill(missing, max_len - min_len_2)),
        C_tw2 = vcat(avg_c_tw2, fill(missing, max_len - min_len_2)),
    )

    # 4. Write to the file
    CSV.write(parsed_args["output_file"], df)

    println("Results saved to $(parsed_args["output_file"])")
end

let
    main()
end

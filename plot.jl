using CSV
using DataFrames
using Plots

# --- Main Plotting Function ---
function create_plot()
    # 1. Define the input file from your simulation
    input_file_name = "test_results_glassy_new"
    input_file = input_file_name * ".csv" # Make sure this matches your output file!
    output_plot_file = input_file_name * ".png"
    println("Loading data from $input_file")

    # 2. Load the CSV file into a DataFrame
    #    The 'try-catch' block handles errors if the file isn't found
    df = try
        CSV.read(input_file, DataFrame)
    catch e
        println("Error: Could not read file '$input_file'.")
        println("Did you run the simulation yet?")
        println(e)
        return # Exit the function
    end

    # 3. Extract the data.
    #    We need to filter out the 'missing' values we added as padding.

    # Get data for the first wait time (t_w1)
    # dropmissing() creates a new, clean DataFrame with only the valid rows
    df_tw1 = dropmissing(df, [:delta_t_1, :C_tw1])

    # Get data for the second wait time (t_w2)
    df_tw2 = dropmissing(df, [:delta_t_2, :C_tw2])

    println("Data loaded. Creating plot...")

    # 4. Create the plot

    # Create the initial plot object
    # We use plot() for the first curve
    p = plot(
        df_tw1.delta_t_1,  # X-axis data
        df_tw1.C_tw1,      # Y-axis data
        label = "t_w1", # Assumes t_w1 is the first time
        xlabel = "Time Difference (t' - t_w)",
        ylabel = "Autocorrelation C(t_w, t')",
        title = output_plot_file,
        legend = :topright,
        linewidth = 2,
    )

    # Add the second curve to the *same* plot
    # We use plot!() (with a '!') to modify the existing plot 'p'
    plot!(
        p,                 # The plot to modify
        df_tw2.delta_t_2,  # X-axis data
        df_tw2.C_tw2,      # Y-axis data
        label = "t_w2", # Assumes t_w2 is the first time
        linewidth = 2,
    )

    # Optional: Use a log scale for the x-axis to see the decay better
    # You can comment this out if you prefer a linear scale.
    #    plot!(p, xscale = :log10)

    # 5. Save the plot to a file
    savefig(p, output_plot_file)

    println("Plot saved to $output_plot_file")
end

# --- Run the plotting function ---
create_plot()

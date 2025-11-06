using CSV
using DataFrames
using Plots
using ArgParse

function parse_plot_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_file", "-i"
            help = "Path to the input CSV file"
            arg_type = String
            default = "results.csv"
        "--output_file", "-o"
            help = "Path to save the output plot"
            arg_type = String
            default = "energy_plot.png"
    end
    return parse_args(s)
end

function create_energy_plot(args)
    input_file = args["input_file"]
    output_file = args["output_file"]

    println("Loading data from $input_file...")
    df = try
        CSV.read(input_file, DataFrame)
    catch e
        println("Error: Could not read file '$input_file'.")
        println(e)
        return
    end

    if isempty(df)
        println("Error: Data file is empty.")
        return
    end

    println("Data loaded. Creating plot...")

    # Plot 1: Energy vs. log(time) (to see plateaus)
    p1 = plot(
        df.time,
        df.avg_energy,
        xlabel = "Time (t)",
        ylabel = "Average Energy E(t)",
        title = "Energy Decay (Log Time Axis)",
        xscale = :log10, # <-- This is the important part
        legend = false,
        linewidth = 2
    )

    # Plot 2: Energy vs. linear(time)
    p2 = plot(
        df.time,
        df.avg_energy,
        xlabel = "Time (t)",
        ylabel = "Average Energy E(t)",
        title = "Energy Decay (Linear Time Axis)",
        legend = false,
        linewidth = 2
    )

    # Combine the two plots into one figure
    p_final = plot(p1, p2, layout = (2, 1), size = (800, 1000))

    savefig(p_final, output_file)
    println("Plot saved to $output_file")
end

# --- Run the main function ---
args = parse_plot_args()
create_energy_plot(args)
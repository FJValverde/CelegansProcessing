#! /usr/local/bin/julia
using DrWatson
@quickactivate "CelegansProcessing"

#%% Import packages
using CSV
using DataFrames

using Revise#This is to be able to modify module Celegans and reload it
using Celegans

#%% Ignored packages for data loading and preprocessing
# using XLSX
# using ColorSchemes
# using Colors
# using DiffEqBase
# using DifferentialEquations
# using LinearAlgebra
# using ModelingToolkit

# FVA: so far 09/08/23

print("1. Loading all information dataframes...")
data_connect_synaptic = CSV.read(datadir("exp_pro","data_connect_synaptic.csv"), DataFrame);
data_connect_gap = CSV.read(datadir("exp_pro","data_connect_gap.csv"), DataFrame);
data_connect_monoamine = CSV.read(datadir("exp_pro","data_connect_monoamine.csv"), DataFrame);
data_connect_neuropeptide = CSV.read(datadir("exp_pro","data_connect_neuropeptide.csv"), DataFrame);
println("Done!")


# 2. EDA: exploratory data analysis
using Plots

using Celegans: nNeurons

#%% Data exploration 
println("2.1. Plot the number of connections of the synaptic and gap coonectomes")


# 2.1.1 Obtaining the synaptic connectome
(synaptic_number, synaptic_connections) =
    connectomeToMultivaluedMatrix(data_connect_synaptic) 
# Check that they do not return void matrices
sum(synaptic_connections)
sum(synaptic_number)

# synaptic_number = zeros(302, 302);      # To store the number of connections
# synaptic_connections = zeros(302, 302); # To store if there is a connection
# # Iterate through the datatframe and add the information to the corresponding point
# for link in eachrow(data_connect_synaptic)
#     from_index = link[:IndexSending]
#     to_index = link[:IndexReceiving]
#     #value = link[:Number]
#     synaptic_number[from_index, to_index] = link[:Number]#value
#     synaptic_connections[from_index, to_index] = 1 #FVA: should this be boolean?
# end
# Plot both matrices
spy(synaptic_number, plot_title= "Number of synaptic connections", xlabel = "Sending neuron index", ylabel = "Rceiving neuron index", color=:berlin)
spy(synaptic_connections, plot_title= "Synaptic connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")


#%% For viewing the information of the frequency of the number of connections
data_s = vec(synaptic_number);           # Flatten the matrix into a 1D array
data_s = filter(x -> x != 0, data_s);    # Take out the values that are 0 for a true result
# Generate histogram
histogram(data_s, bins=:scott, 
    xticks = unique(data_s), legend = false, normalize=:probability,
          xlabel = "Number of connections between two neurons",
          ylabel = "Frequency", title = "Synaptic frequency histogram")
histogram!(size=(760,500))


println("2.1.2 Obtaining the gap junction connectome")
# Repeat the previous process for GAP junctions
# gap_number = zeros(302, 302)
# gap_connections = zeros(302,302)
# for link in eachrow(data_connect_gap)
#     from_index = link[:IndexSending]
#     to_index = link[:IndexReceiving]
#     value = link[:Number]
#     gap_number[from_index, to_index] = value
#     gap_connections[from_index, to_index] = 1
# end
(gap_number, gap_connections) = connectomeToMultivaluedMatrix(data_connect_gap)

# Exploration
spy(gap_number, plot_title= "Number of gap connections among neurons",
    xlabel = "Sending neuron index",
    ylabel = "Receiving neuron index",
    color=:berlin)
spy(gap_connections, plot_title= "Gap connections among neurons",
    xlabel = "Sending neuron index",
    ylabel = "Receiving neuron index")
# FVA.HYPOTHESIS: the GJ connectome is symmetrical.
if (gap_connections != transpose(gap_connections))
    println("WARNING: The gap junction connectome is not symmetrical!")
end

# Histogram for the number of gap junctions between neurons.
data = vec(gap_number);
data = filter(x -> x != 0, data);
# Generate histogram
histogram(data, bins=:scott, 
    xticks = unique(data), legend = false, normalize=:probability,
    xlabel = "Number of connections between two neurons", ylabel = "Frequency", title = "Gap frecuency histogram")
histogram!(size=(800,500))


print("ESTO ES EL SURPRISAL (NO ESTA BIEN TODAVIA)...")
# Calculate surprisal values for the data
surprisal_data = -log2.(data);

# Sort the data and surprisal values
sorted_data, sorted_surprisal = sort(data), surprisal_data[sortperm(data)]

# Plot the surprisal curve
plot!(sorted_data, sorted_surprisal, 
    linecolor=:blue, linewidth=2, 
    xlabel="Values", ylabel="Surprisal", title="Surprisal Curve")
#gap_number = log.(gap_number .+ 0.7)
println("AQUI TERMINA EL SURPRISAL")


println("2.1.3 Obtaining the monoamine connectome")
#println("Plot the connections by Monoamines")
# Create a matrix for each type
# mono_connections_tyr = zeros(nNeurons,nNeurons);
# mono_connections_oct = zeros(nNeurons,nNeurons);
# mono_connections_dop = zeros(nNeurons,nNeurons);
# mono_connections_ser = zeros(nNeurons,nNeurons);
# # Add each connection to the corresponding matrix
# for link in eachrow(data_connect_monoamine)
#     from_index = link[:IndexSending]
#     to_index = link[:IndexReceiving]
#     if link[:Type] == "tyramine"
#         mono_connections_tyr[from_index, to_index] = 1
#     elseif link[:Type] == "octopamine"
#         mono_connections_oct[from_index, to_index] = 1
#     elseif link[:Type] == "dopamine"
#         mono_connections_dop[from_index, to_index] = 1
#     elseif link[:Type] == "serotonin"
#         mono_connections_ser[from_index, to_index] = 1
#     end
# end

mono_connections_tyr =
    connectomeToBinaryMatrix(
        filter(:Type => ==("tyramine"),data_connect_monoamine)
    );
mono_connections_oct = 
    connectomeToBinaryMatrix(
        filter(:Type => ==("octopamine"),data_connect_monoamine)
    );
mono_connections_dop =
    connectomeToBinaryMatrix(
        filter(:Type => ==("dopamine"),data_connect_monoamine)
    );
mono_connections_ser =
    connectomeToBinaryMatrix(
        filter(:Type => ==("serotonin"),data_connect_monoamine)
    );

# Plot all in the same graph
" FVA: these plots have to be extended in the y range so that all
vertical neuron indices appear."
spy(mono_connections_tyr, color = :orange,
    plot_title= "Monoamine connections among neurons",
    xlabel = "Sending neuron index", ylabel = "Receiving neuron index")
spy!(mono_connections_oct, color = :red)
spy!(mono_connections_dop, color = :green)
spy!(mono_connections_ser, color = :blue,legend=:outertopright)


println("2.1.4 Plot the connections by Neuropeptides (FALTA SABER QUE SIGNIFICA CADA COLUMNA)")
# Create one matrix for each column of the dataframe
# neuropep_connections_Type1 = zeros(302,302)
# neuropep_connections_Type2 = zeros(302,302)
# # Append the value to the matrix
# for link in eachrow(data_connect_neuropep)
#     from_index = link[:IndexSending]
#     to_index = link[:IndexReceiving]
#     if link[:Type1] !== nothing
#         neuropep_connections_Type1[from_index, to_index] = 1
#     end
#     if link[:Type2] !== nothing
#         neuropep_connections_Type2[from_index, to_index] = 1
#     end
# end
neuropep_connections_Type1 =
    connectomeToBinaryMatrix(
        filter(:Type1 => !=(nothing), data_connect_neuropeptide)
    );
neuropep_connections_Type2 =
    connectomeToBinaryMatrix(
        filter(:Type2 => !=(nothing), data_connect_neuropeptide)
    );

spy(neuropep_connections_Type1, color = :blue,
    plot_title= "Neuropeptides connections among neurons",
    xlabel = "Neuron index", ylabel = "Neuron index")
spy!(neuropep_connections_Type2, color = :red)

# FVA: We should distinguish the type of neuroamine in this map.
# FVA: Q. Which of these plots should be saved for explanation?
print("3. Saving the matrix data for later modelling...")

println("Done!")

println("4. Environment description")
using Pkg;Pkg.status()

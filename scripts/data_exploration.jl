#! /usr/local/bin/julia
#
# Author: FVA, based on V. Bokach's long monolithic script.
# 
# This script has two intended usages:
# - as an interactive script, to visualise the multiple layers of the
# connectome of C.elegans, and generate the matrices related to those. 
# - as a batch script to generate said matrices.
#
# TODO: transform this into a jupyter or Quarto document.
using DrWatson# Scientific project management with Julia

@quickactivate "CelegansProcessing"

#%% Import packages
using CSV
using DataFrames
using TidierData
using TidierStrings
using Plots
using JLD2
using SparseArrays

using Revise
#This is to be able to modify module Celegans and reload it.
# CAVEAT: use it only during Celegans development.
using Celegans

# 2. EDA: exploratory data analysis
# 2.1 A forward dictionary from names of neurons to indices


# FVA: the next is better for data exploration?
#%% Data homogeneization
print("TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR GAP AND SYNAPTIC CONNECTIONS...")

# Solution by FVA:

"""

        neuron_list

A dataframe describing what we know about individual C.elegans neurons. Generated by the preprocesing script and re-read from Celegans.Files.neuronList

"""
neuron_list =
    @chain CSV.read(datadir("exp_pro",Celegans.Files.neuronList),
                    DataFrame;
                    stringtype=String) begin
      @unite(NT, [Neurotransmitter1, Neurotransmitter2], " ")
      @mutate(NType = case_when(
                        str_detect(NT, "GABA") => "Inhibitory",
                        true  => "Excitatory"))
    end;
describe(neuron_list)
@assert !all(neuron_list.NType .== "Excitatory")
@assert length(Set(neuron_list.Neuron)) == Celegans.nNeurons
@assert reduce(max, neuron_list.Index) == Celegans.nNeurons
# There are both inhibitory and excitatory neurons in Celegans. 
@assert Set(unique(neuron_list.NType)) == Set(["Inhibitory", "Excitatory"])
@assert sum(neuron_list.NType .== "Inhibitory") != 0

# a) Create dictionary for neuron names to sorted numbers
indexByName = Dict(
    zip(neuron_list.Neuron,neuron_list.Index)
);

"""
        getIndexByName

A function to obtain the index of a neuron in a standard ranking by distance to the tip. It checks that the name exists or will produce 'nothin'

Use as indexByName["BAGL"], but checks that the name is bound.
@example

indexByName["BAGL"] == getIndexByName("BAGL")
isnothing(getIndexByName("BAGL")) 

"""
function getIndexByName(NeuronName::String)
    return(get(indexByName,NeuronName,nothing))
end

"""

        nameByIndex

An inverse index-to-names for each neuron. Maintained as an array and thus accessible.
"""
nameByIndex = neuron_list.Neuron;
#These two are inverses. (mapreduce is more elegant)
@assert all(map(x -> nameByIndex[indexByName[x]] == x,nameByIndex))
@assert all(map(x -> indexByName[nameByIndex[x]] == x,1:nNeurons))

#FVA.TODO: THis should go in the initialization of the module Celegans.

print("1. Loading all information dataframes and transforming neuron names...")
data_connect =
    @chain CSV.read(datadir("exp_pro","data_connect.csv"),
                    DataFrame; stringtype=String) begin
        @filter(!(Type=="NMJ"))
        @mutate(IndexSending = getIndexByName(Sending),
                IndexReceiving = getIndexByName(Receiving))
    end;
println(describe(data_connect))
println("Count of Gap Junctions and Chemical Synapses")
data_connect_grouped = 
    @chain data_connect begin
        #arrange(Type)
        groupby(:Type, sort=true)#Groups and sorts at the same time
    end;
@summarize(data_connect_grouped, Count = sum(Number))


#data_connect_gap = CSV.read(datadir("exp_pro","data_connect_gap.csv"), DataFrame);
data_connect_monoamine =
    @chain CSV.read(datadir("exp_pro","data_connect_monoamine.csv"),
                    DataFrame; stringtype=String) begin
        @mutate(IndexSending = getIndexByName(Sending),
                IndexReceiving = getIndexByName(Receiving))
    end;
println(describe(data_connect_monoamine))
println("Counts for MONOAMINE mediated connections...")
@chain data_connect_monoamine begin
    @group_by(Neurotransmitter)
    @summarize(Count=n())
end


data_connect_neuropeptide =
    @chain CSV.read(datadir("exp_pro","data_connect_neuropeptide.csv"),
                    DataFrame; stringtype=String) begin
        @mutate(IndexSending = getIndexByName(Sending),
                IndexReceiving = getIndexByName(Receiving))
    end;
println(describe(data_connect_neuropeptide))
@chain data_connect_neuropeptide begin
    @group_by(Neurotransmitter)
    @summarize(Count=n())
end
println("Done!")


#%% Data exploration 
println("The types of connetions are:", unique(data_connect.Type))
describe(data_connect)
# 2.1.1 Obtaining the GAP JUNCTION CONNECTOME
data_connect_gap = data_connect_grouped[(Type="EJ",)];#
@assert unique(data_connect_gap.Type) == ["EJ"]
data_connect_gap =
    @chain data_connect_gap begin
        @group_by(IndexSending,IndexReceiving)
        #@summarise(Count=sum(Number))
        combine(:Number => sum,
                renamecols=false,#Do not rename :Number
                #copycols=true, #all information extant required
                ungroup=true) #But return a non-grouped DF.
    end;
@assert sum(data_connect_gap.Number) == sum(data_connect_grouped[(Type="EJ",)].Number)
println("Number of gap junctions: $(sum(data_connect_gap.Number))")
gap_connectome = sparse(
    data_connect_gap.IndexReceiving,#Post synaptic neuron
    data_connect_gap.IndexSending,#Pre-synaptic neuron
    data_connect_gap.Number,#Number of neurons described
    nNeurons,nNeurons)
@assert sum(gap_connectome) == sum(data_connect_grouped[(Type="EJ",)].Number) "Warning: modelled number of GJ is different from recorded data!"
# FVA.HYPOTHESIS: the GJ connectome is symmetrical.
# 1. Find asymmetrical component using standard theory
asym_gap_connectome = (gap_connectome - gap_connectome')/2;
# 1.1. Visualize it
spy(asym_gap_connectome,
    plot_title= "Asymmetric part of GJ connectome",
    ylabel = "Sending neuron index",
    xlabel = "Receiving neuron index",
    color=:berlin)
# 1.1. How many entries are non-symmetrical and what are their values?
if  nnz(asym_gap_connectome) != 0
    # A.FVA: here we symmetrize it:
    print("The gap connectome is not symmetrical in $(nnz(asym_gap_connectome)) values:")
    println(findnz(asym_gap_connectome))
    # 2. Symmetrize the gap junction connectome
    print("About to symmetrize the GJ connectome...")
    gap_connectome = gap_connectome - asym_gap_connectome;
    println("Done!")
    # else
    # Exploration: what happens if the GJ connectome is not symmetrical?
    # end
end
@assert all(gap_connectome == gap_connectome') "WARNING: The gap junction connectome is not symmetrical!"

println("2.1.2 Obtaining the synaptic connectome ")
"""
data_connect_synaptic

Select only synaptic connection data and add over all possible combinations of
(Sending, Receiving) columns. 
"""
data_connect_synaptic = data_connect_grouped[(Type="S",)]
@assert unique(data_connect_synaptic.Type) == ["S"]# "EJ" < "S"
# 0. Loaded and aggregated by key (pre-, post-neuron)
data_connect_synaptic = 
    @chain data_connect_synaptic begin
        @group_by(IndexSending,IndexReceiving)
        #@summarise(Count=sum(Number))
        combine(:Number => sum,
                renamecols=false,#Do not rename :Number
                #copycols=true, #all information extant required
                ungroup=true) #But return a non-grouped DF.
    end;
describe(data_connect_synaptic)
@assert sum(data_connect_synaptic.Number) == sum(data_connect_grouped[(Type="S",)].Number)
#@filter(data_connect, Type == "S")
println("Number of chemical synapses: $(sum(data_connect_synaptic.Number))")

# 1. The synaptic data are enriched with the type of the synapse
data_connect_synaptic_all =
    @chain data_connect_synaptic begin
        leftjoin(
            select(neuron_list,
                   [:Index,
                    #:Neurotransmitter1,
                    #:Neurotransmitter2,
                    :NType]),
            on = :IndexSending => :Index
        )#keep synapses even if no info
        #@select(-(Neurotransmitter1:Neurotransmitter2))
        @group_by(NType)#Prepare split
    end;
@chain data_connect_synaptic_all begin
    @ungroup()
    describe()
end
@chain data_connect_synaptic_all begin
    @summarise(Count=sum(Number))
end
#TODO: Check that these number make sense in the literature

# 2. Split by inhibitory and excitatory parts
data_connect_synaptic_inh = 
    data_connect_synaptic_all[(NType = "Inhibitory",)]
# The complement (yes, as naive as this) are the excitatory
data_connect_synaptic_exc  =
    data_connect_synaptic_all[(NType = "Excitatory",)]
#@assert sum(data_connect_synaptic_excitatory.Neurotransmmitter1 == "GABA") == 0 

# # 2. Add a column declaring each neuron of :NType Inhibitory or Excitatory
# Set(unique(neuron_list.Neurotransmitter1))
# sum(neuron_list.Neurotransmitter1 .== "GABA") + 
#     sum(skipmissing(neuron_list.Neurotransmitter2) .== "GABA")
# sum(neuron_list.Neurotransmitter1 .== "Unknown")
# sum(neuron_list.Neurotransmitter1 .== "Monoamine")
# # 3. Group on type of neuron

# 3. For each group, generate a connectome
#As per Kunert et als, 
inhibitory_connectome = sparse(
    data_connect_synaptic_inh.IndexReceiving,#Post synaptic neuron
    data_connect_synaptic_inh.IndexSending,#Pre-synaptic neuron
    data_connect_synaptic_inh.Number,#Number of neurons described
    nNeurons,nNeurons)
excitatory_connectome = sparse(
    data_connect_synaptic_exc.IndexReceiving,#Post synaptic neuron
    data_connect_synaptic_exc.IndexSending,#Pre-synaptic neuron
    data_connect_synaptic_exc.Number,#Number of neurons described
    nNeurons,nNeurons)
#Assertion: synaptic connections are distributed among inhibitory and excitatory connections. 
@assert nnz(inhibitory_connectome) + nnz(excitatory_connectome) == nrow(data_connect_synaptic)
# Assertion: total synaptic strength is distributed between inhibitory and excitatory connections.
@assert sum(inhibitory_connectome) + sum(excitatory_connectome) == sum(data_connect_synaptic.Number)

# # Looking for repeated pairs
# @assert nnz(synaptic_number) ==
#     length(unique(zip(data_connect_synaptic.IndexSending,
#                   data_connect_synaptic.IndexReceiving)))
# #In conclusions, there are repeated entries: are they from names?
# @assert nnz(synaptic_number) ==
#     length(unique(zip(data_connect_synaptic.Sending,
#                       data_connect_synaptic.Receiving)))
# #YES! Originally the pairs of neurons are not mutually exclusive!


# Exploring the synaptic connectome
gr()# Invoke the primitives

#%% For viewing the information of the frequency of the number of connections
# a) Exploring the EXCITATORY connectome
# Hypothesis: It's a decreasing exponential
# TODO: estimate the constant.
# data_s = vec(excitatory_connectome);# Flatten the matrix into a 1D array
# data_s = filter(x -> x != 0, data_s);
# # TODO: use chain notation on data processing here too!
data_s = 
    @chain vec(excitatory_connectome) begin# Flatten the matrix into a 1D array
        # Take out the values that are 0 for a true result, because they dominate the histogram!
        filter(x -> x != 0, _)
        #map(log,_)
    end;
# Recall that this is a distribution over counts (a multinouilli?).
# Generate histogram
histogram(data_s,
          bins=:scott, 
          xticks = unique(data_s), legend = false, normalize=:probability,
          xlabel = "Number of connections between two neurons",
          ylabel = "Frequency",
          title = "Synaptic frequency histogram")
histogram!(size=(760,500))

# Plot both matrices
spy(excitatory_connectome,
    plot_title= "Number of excitatory synaptic connections",
    xlabel = "Receiving neuron index",
    ylabel = "Sending neuron index",
    #color=:berlin)#too bland in difference between H and L.
    #color=:thermal)#Highest colour: yellow, too similar to white.
    #color=:bluesreds)#too concentrated in low values... Let's use log
    color=cgrad(:darkrainbow, scale = :log))#, rev=true))#
spy(inhibitory_connectome,
    plot_title= "Number of inhibitory synaptic connections",
    xlabel = "Receiving neuron index",
    ylabel = "Sending neuron index",
    #color=:berlin)
    color=cgrad(:darkrainbow, scale = :log))#, rev=true))#
#spy(synaptic_connections, plot_title= "Synaptic connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")

println("2.1.2 Explring the gap junction connectome")
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
#(gap_number, gap_connections) = connectomeToMultivaluedMatrix(data_connect_gap)

# Histogram for the number of gap junctions between neurons.
# Conclusion: looks like a decreasing exponential.
# TODO: infer the constant
# data = vec(gap_number);
# data = filter(x -> x != 0, data);
data_gj =
    @chain vec(gap_connectome) begin
        filter(x -> x != 0, _)
    end;
# Generate histogram
histogram(data_gj,
          bins=:scott, 
          xticks = unique(data_gj), legend = false, normalize=:probability,
          xlabel = "Number of connections between two neurons",
          ylabel = "Frequency",
          title = "Gap Junction frequency histogram")
histogram!(size=(800,500))

println("Exploring the fill pattern of the GJ connectome...")
spy(gap_connectome,
    plot_title= "Number of gap connections among neurons",
    ylabel = "Sending neuron index",
    xlabel = "Receiving neuron index",
    color=:berlin)
# spy(gap_connections, plot_title= "Gap connections among neurons",
#     xlabel = "Sending neuron index",
#     ylabel = "Receiving neuron index")
println("Done!")

# # print("ESTO ES EL SURPRISAL (NO ESTA BIEN TODAVIA)...")
# # Calculate surprisal values for the data
# surprisal_data = -log2.(data);

# # Sort the data and surprisal values
# sorted_data, sorted_surprisal = sort(data), surprisal_data[sortperm(data)]

# # Plot the surprisal curve
# plot!(sorted_data, sorted_surprisal, 
#     linecolor=:blue, linewidth=2, 
#     xlabel="Values", ylabel="Surprisal", title="Surprisal Curve")
# #gap_number = log.(gap_number .+ 0.7)
# println("AQUI TERMINA EL SURPRISAL")


#println("2.1.3 Obtaining the monoamine connectome")
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


# mono_connections_tyr =
#     connectomeToBinaryMatrix(
#         filter(:Type => ==("tyramine"),data_connect_monoamine)
#     );
# mono_connections_oct = 
#     connectomeToBinaryMatrix(
#         filter(:Type => ==("octopamine"),data_connect_monoamine)
#     );
# mono_connections_dop =
#     connectomeToBinaryMatrix(
#         filter(:Type => ==("dopamine"),data_connect_monoamine)
#     );
# mono_connections_ser =
#     connectomeToBinaryMatrix(
#         filter(:Type => ==("serotonin"),data_connect_monoamine)
#     );

# # Plot all in the same graph
# " FVA: these plots have to be extended in the y range so that all
# vertical neuron indices appear."
# spy(mono_connections_tyr, color = :orange,
#     plot_title= "Monoamine connections among neurons",
#     xlabel = "Sending neuron index",
#     ylabel = "Receiving neuron index")
# spy!(mono_connections_oct, color = :red)
# spy!(mono_connections_dop, color = :green)
# spy!(mono_connections_ser, color = :blue,legend=:outertopright)


# println("2.1.4 Plot the connections by Neuropeptides (FALTA SABER QUE SIGNIFICA CADA COLUMNA)")
# # Create one matrix for each column of the dataframe
# # neuropep_connections_Type1 = zeros(302,302)
# # neuropep_connections_Type2 = zeros(302,302)
# # # Append the value to the matrix
# # for link in eachrow(data_connect_neuropep)
# #     from_index = link[:IndexSending]
# #     to_index = link[:IndexReceiving]
# #     if link[:Type1] !== nothing
# #         neuropep_connections_Type1[from_index, to_index] = 1
# #     end
# #     if link[:Type2] !== nothing
# #         neuropep_connections_Type2[from_index, to_index] = 1
# #     end
# # end
# neuropep_connections_Type1 =
#     connectomeToBinaryMatrix(
#         filter(:Type1 => !=(nothing), data_connect_neuropeptide)
#     );
# neuropep_connections_Type2 =
#     connectomeToBinaryMatrix(
#         filter(:Type2 => !=(nothing), data_connect_neuropeptide)
#     );

# spy(neuropep_connections_Type1, color = :blue,
#     plot_title= "Neuropeptides connections among neurons",
#     xlabel = "Neuron index", ylabel = "Neuron index")
# spy!(neuropep_connections_Type2, color = :red)

# FVA: We should distinguish the type of neuroamine in this map.
# FVA: Q. Which of these plots should be saved for explanation?

print("3. Saving the matrix data for later modelling...")
#@write datadir(datadir("exp_pro","connectome_matrices.hdf5")) data_connect_gap data_connect_synaptic data connect_monoamine data connect_neuropeptide
# Despite claims to the contrary, HDF5 does not define @safe @load as Matlab does.
#using JLD2
#cfile = datadir("exp_pro",Celegans.Files.MlConnectome);
#Export this name to a submodule of Celegans.Paths.ToMLConnectome

jldsave(datadir("exp_pro",Celegans.Files.mlConnectome);
        nameByIndex = nameByIndex,#load with global nameByIndex...
        indexByName = indexByName,#load with global nameByIndex...
        neuronList = neuron_list,#Extended neuron list with type of neuron
        dataGapJunctions=data_connect_gap,#DataFrame
        gapJunctionConnectome=gap_connectome,#Sparse matrix

        dataSynapticAll=data_connect_synaptic_all,#DataFrame
        InhibitoryConnectome=inhibitory_connectome,#Sparse matrix
        ExcitatoryConnectome=excitatory_connectome,#Sparse matrix

        dataMonoamine=data_connect_monoamine,#DataFrame
        # monoTyramineBinary=mono_connections_tyr,
        # monoOctopamineBinary=mono_connections_oct,
        # monoDopamineBinary=mono_connections_dop,
        # monoSerotonineBinary=mono_connections_ser,

        dataNeuropeptide=data_connect_neuropeptide#DataFrame
        )
println("Done!")
#load with:
# connectomes = load(datadir("exp_pro",Celegans.Files.mlConnectome))
# data_connect_gap = connectomes["data_connect_gap"]
#or load(cfile; data_connect_gap)

# TODO: see how to properly export nested modules in Julia
jldsave(datadir("exp_pro",Celegans.Files.dictionaries);
        indexByName)
#TODO: store the dictionary of neurotransmitters.

println("4. Environment description")
using Pkg;Pkg.status()


#FVA SO FAR DEBUGGING HERE! 
exit(0)
# SO FAR 29/08/23

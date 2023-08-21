#! /usr/local/bin/julia
using DrWatson
@quickactivate "CelegansProcessing"

#%% Import packages
using CSV
using DataFrames
using TidierData
using TidierStrings
using XLSX
#using JLD2

#%% Ignored packages for data loading and preprocessing
# using ColorSchemes
# using Colors
# using DiffEqBase
# using DifferentialEquations
# using LinearAlgebra
# using ModelingToolkit

using Revise
using Celegans


#%% DATA IMPORT
print("0. DATA IMPORT: ")
println("extracting the data from different files and stores as dataframes...")

print("* DATA for pharingeal neuron connectivity...")
"""
Pharingeal neuron connectivity data

Sending: Name of sending neuron
Receiving: Name of receiving neuron
Number: Number of synapses between the given neuron pair.
Type: Type of synapse: S: synaptic; G: gap.
"""
data_connect_phar =
    @chain CSV.read(datadir("exp_raw","ConexionsPharyngeal.csv"),
                    DataFrame;
                    stringtype=String) begin
        #select([:Sending,:Receiving,:Type,:Number])#reordercols
        @rename(OldType = Type)
        @mutate(Type = if_else(OldType == "G", "EJ", OldType))#shared codes
        #@select(Sending, Receiving, Type, Number)
        @select(-(OldType))#dispose of dummy
    end
# data_connect_phar = CSV.read("RawData/ConexionsPharyngeal.csv", Data
# Frame)
describe(data_connect_phar)
println("Done!")
@assert Set(unique(data_connect_phar[!,:Type])) == Set(["S", "EJ"])

print("* DATA for the connectome of other neurons...")
"""
These data come from WormAtlas in NeuronConnect.xls(x)

N1: Neuron 1 name
N2: Neuron 2 name
Type: Type of synapse: 
    S: Send or output (Neuron 1 pre-synaptic to Neuron 2); 
    Sp: Send-poly (Neuron 1 is pre-synaptic to more than one postsynaptic partner. Neuron 2 is just one of these post-synaptic neurons); 
    R: Receive or input (Neuron 1 is post-synaptic to Neuron 2); 
    Rp: Receive-poly (Neuron 1 is one of several post-synaptic partners of Neuron 2.); 
    EJ: Electric junction; 
    NMJ: Neuromuscular junction (only reconstructed NMJ's are represented).
It must be noted that at polyadic synaptic sites, not all “send-poly” were faithfully labeled 
    as such in White et al, 1986. Some pre-synaptic connections were labeled simply as “sends”. 
    Reconciliation of chemical synapses did not previously distinguish between send from send-poly 
    and receive from receive-poly. In this new reconciliation, the total number of send and send-poly 
    is equal to the total number of receive and receive-poly (S+Sp=R+Rp). Every documented synapse is 
    now listed in this Table, both with respect to the sending neuron and with respect to the receiving neuron(s).
Nbr: Number of synapses between the given neuron pair.
"""
data_connect_neuron =
    @chain DataFrame(XLSX.readtable(
        datadir("exp_raw","NeuronConnect.xlsx"),"Sheet1"; infer_eltypes=true)) begin
        rename([:Sending,:Receiving,:OldType,:Number])#Standard names
        @mutate(Type = if_else(OldType == "Sp", "S", OldType))#FVA: Send multiple is also a Send.
        @select(-(OldType))#Delete original type
    end;
# data_connect_neuron = DataFrame(XLSX.readtable("RawData/NeuronConnect.xlsx", "Sheet1"))
describe(data_connect_neuron)
println("Done!")
@assert Set(unique(data_connect_neuron[!,:Type])) == Set(["S", "EJ", "NMJ", "R","Rp"])


print("* DATA for the connectome of ALL neurons...")
data_connect =
    @chain vcat(data_connect_phar, data_connect_neuron, 
                cols = [:Sending,:Receiving,:Type, :Number]) begin
        @filter(!(Type in ("R","Rp") ))#keep only S, EJ, NMJ
    end;
println(describe(data_connect))
println("Done!")
@assert Set(unique(data_connect[!,:Type])) == Set(["S", "EJ", "NMJ"])

print("1.2 READ CONNECTIONS FOR MONOAMINES...")
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type: Name of the monoamine
Monoamine: MA stands for monoamine
Specific: The specific type of monoamine
"""
data_connect_monoamine =
    @chain CSV.read(datadir("exp_raw","MonoaminesConnect.csv"), DataFrame;
             types=String) begin
        @rename(Sending = Neuron1, Receiving = Neuron2, Neurotransmitter = Type)
        @select(-(Monoamine))
        @mutate(Type="MA")#For "MonoAmine"         
    end
# data_connect_monoamine = CSV.read("RawData/MonoaminesConnect.csv", DataFrame)
println(describe(data_connect_monoamine))
println("Done!")


# DEBUGGING: SO far 21/08/23
print("1.3 READ THE CONNECTIONS FOR neuropeptides...")
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type1: One neuropeptide
Type2: Another neuropeptide
"""
data_connect_neuropep =
    @chain CSV.read(datadir("exp_raw","NeuropeptidesConnect.csv"), DataFrame;types=String) begin
        @rename(Sending = Neuron1, Receiving = Neuron2)
        @mutate(Type = "NP", Neurotransmitter = string(Type1, "|", Type2))
        @select(-(Type1))#, Type2))
        @select(-(Type2))#FVA: funny that this cannot be made for multiple columns
    end;
# data_connect_neuropep = CSV.read("RawData/NeuropeptidesConnect.csv", DataFrame)
println(describe(data_connect_neuropep))
println("Done!")


print("1.0 Ordered list of neurons sorted by distance to tip, giving a standard rank on them. ")
"""
        Data on individual neurons

These data come from WormAtlas in NeuronType.xls(x)(sheet1)

Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
Span: Length of neuron span. Neurons spanning <25% of worm body (e.g., motor neurons in the ventral cord, 
    neurons with processes confined to the nerve ring and neurons confined in the mid-body) are defined to 
    have short spans (S). All other neurons are defined to have long spans (L).
Ambiguity: If applicable, code for the type of ambiguity. Codes beginning with M denote ambiguity 
    citedin White et al, 1986. Codes beginning with R denote ambiguity found in reconstructions during 
    update of wiring diagram (MB=cell body position ambiguous, MTS=tail synapses ambiguous and/or sparse 
    connections in the tail; MAS=anterior body ambiguous and/or sparse connections in the anterior; MD=dorsal 
    side ambiguous; MAD=anterior and dorsal side ambiguous; MS=neurons with sublateral processes not covered 
    by reconstructions. RDI=dorsal reconstruction incomplete; RDM=dorsal reconstruction completely missing; 
    RVI=ventral reconstruction incomplete.)
TotHead: Total number of synapses in the head including EJ and NMJ.
TotTail: Total number of synapses in the tail including EJ and NMJ.
TotMid: Total number of synapses in the mid-body including EJ and NMJ.
S_Head: Number of “sends” or output synapses in the head, includes send polyadic synapses (see Figure 1) and NMJ.
R_Head: Number of “receives” or input synapses (includes polyadic synapses) in the head.
S_Mid: Number of “sends” or output synapses in the mid-body, includes polyadic synapses and NMJ.
R_Mid: Number of “receives” or input synapses (includes polyadic synapses) in the mid-body.
S_Tail: Number of “sends” or output synapses in the tail, includes polyadic synapses and NMJ.
R_Tail: Number of “receives” or input synapses (includes polyadic synapses) in the tail.
AY NeuronType: Letter codes denoting ganglion group as defined by Achacoso and Yamamoto W.S., 1991, where 
    A=anterior ganglion, B=dorsal ganglion, C=lateral ganglion, D=ventral ganglion, E=retrovesicular ganglion, 
    F=posterolateral ganglion, G=ventral cord neuron group, H=pre-anal ganglion, J=dorsorectal ganglion, 
    K=lumbar ganglion.
AYNbr: Numeric identifier given by AY for each neuron.
Note:  Sum of S_Head and R_Head does not include electrical junctions (EJ), therefore, does not equal TotHead.  Similar is true for mid-body and tail.
"""
neuron_type =
    DataFrame(XLSX.readtable(datadir("exp_raw","NeuronType.xlsx"),"Sheet1";
                             infer_eltypes = true));
# neuron_type = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet1"))
describe(neuron_type)

"""
        Data on pharingeal neurons

Q.FVA: where do these data come from?

Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
"""
data_type_phar = DataFrame(
    XLSX.readtable(datadir("exp_raw","NeuronType.xlsx"), "Sheet2";
                   infer_eltypes=true)
);
# data_type_phar = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet2"))
describe(data_type_phar)

# FVA: TODO. Maybe this list as well as that of neurotransmitters
# should be stored in the Celegans module, as a data structure. 
#print("1.0 Ordered list of neurons sorted by distance to tip, giving a standard rank on them. ")
#data_neuron_pos = vcat(data_type_phar, neuron_type, cols=:intersect)
unlinked_neuron_pos = # These were provided by VB after some research
    DataFrame(Neuron=["CANL", "CANR"], SomaPosition = 0.61,SomaRegion = "M");
neuron_position_sorted =
    @chain vcat(data_type_phar, neuron_type, cols=:intersect) begin
        rename!([:Neuron,:SomaPosition,:SomaRegion]) 
        @bind_rows(unlinked_neuron_pos) 
    end;
println(describe(neuron_position_sorted))
println("Done!")

print("1.5. Table of neurotransmitters and deduced inhibitory/excitatory activity...")
"""

        Data on neurotransmitter production in neurons

These tables originate from:
- a set of tables found in Pereira et al., 2015  (http://elifesciences.org/content/4/e12432v2),
- combined with information compiled in Loer & Rand, 2016 & 2022, WormAtlas, and
- from Gendrel et al., 2016 and Serrano-Saiz et al., 2017.						

Neuron: Name of the neuron
Neurotransmitter1: Main neurotransmitter
Neurotransmitter2: If the neuron uses another neurotransmitter it is stated here
"""
neuron_by_neurotransmitter =
    @chain DataFrame(
        XLSX.readtable(datadir("exp_raw","Neurotransmitters.xlsx"), "Sheet1";
                       infer_eltypes=true)) begin
    end
# neurotransmitter = DataFrame(XLSX.readtable("RawData/Neurotransmitters.xlsx", "Sheet1"))
describe(neuron_by_neurotransmitter)

# We next return a single table on neurons with their positions and
# neurotransmitters as the natural join of both.
@assert Set(unique(neuron_position_sorted.Neuron)) == Set(unique(neuron_by_neurotransmitter.Neuron))

neuron_list =
    @chain @inner_join(neuron_position_sorted,neuron_by_neurotransmitter,Neuron) begin
        @arrange(SomaPosition)
        @mutate(Index=1:302)#FVA: this fixates the index
    end

                           
# """
# SUMMARY IMPORTING DATA: IMPORTANT VARIABLES: the DataFrames read in...
# 1. data_connect_phar
# 2. data_connect_neuron
# 3. neuron_type
# 4. data_type_phar
# 5. data_connect_monoamine
# 6. data_connect_neuropep
# 7. neurotransmitter

# ...DO NOT TRANSMIT TO OTHER SCRIPTS!
# """
println("Done!")


#%% DATA FRAMES AND CLEANING
println("1. DATAFRAME CLEANING...Done!")


# # VB: From all the dataframe of neuron_type only select the first three columns
# # VB: neuron_type = neuron_type[:, 1:3]
# # 1. Concatenate the two dataframes keeping: :Neuron(name), "Soma Position", "Soma Region"
# # VB: neuron_pos = vcat(neuron_pos_phar, neuron_pos_neuron)
# neuron_pos = vcat(data_type_phar, neuron_type, cols=:intersect)
# # 2. Add the two neurons that do not have any connections
# push!(neuron_pos,["CANL", 0.61, "M"])
# push!(neuron_pos,["CANR", 0.61, "M"])
# # Q.FVA: what is the justification for the soma position of these two added neurons?

# # 4. Sort the dataframe by the position of the neuron in the soma
# neuron_position_sorted = sort!(neuron_pos, [:"Soma Position"])

# # 5. Add a column to know for the future the number of the neuron
# neuron_position_sorted.Index = 1:302
# # FVA: this index based on topological position should be saved, somehow. But also, its inverse!
# println("""
#  IMPORTANT VARIABLES:
# neuron_position_sorted #FVA: As the variable gathering all GJ and Synapsis 
# """)

# 1.1 GAP AND SYNAPTIC
#
# print("1.1 Main GJ and synaptic connectome...")
# FVA: This is massaging the connectome of the pharingeal neurons
# Interchange the "Number" and "Type" columns
#names(data_connect_phar)
#names(data_connect_neuron)
# select!(data_connect_phar, [:Sending,:Receiving,:Type,:Number])
# #FVA: TODO, write the previous select as a TidierData step, like below...
# new_data_connect_phar =
#     @chain data_connect_phar begin
#         select(data_connect_phar, [:Sending,:Receiving,:Type,:Number])
#         @rename(OldType = Type)
#         @mutate(Type = if_else(OldType == "G", "EJ", OldType))
#         @select(Sending, Receiving, Type, Number)
#     end;
# Rename columns so they match the names of "data_connect_phar"
#rename!(data_connect_neuron,[:Sending,:Receiving,:Type, :Number])
# Concatenate the two dataframes about neuron connectivity
# new_data_connect_neuron = 
#     @chain rename(data_connect_neuron,[:Sending,:Receiving,:OldType,:Number]) begin
#         @mutate(Type = if_else(OldType == "Sp", "S", OldType))#FVA: Send multiple is also a Send.
#         @select(Sending, Receiving, Type, Number)
#     end;
# data_connect =
#     @chain vcat(new_data_connect_phar,
#                 new_data_connect_neuron, 
#                 cols = [:Sending,:Receiving,:Type, :Number]) begin
#         @filter(!(Type in ("R","Rp") ))#FVA: keep S, EJ, NMJ
#     end;
# Eliminate the rows in which the "Type" is either receive, receive-poly or NMJ.
#data_connect = data_connect[data_connect.Type .!= "R", :]
#data_connect = data_connect[data_connect.Type .!= "Rp", :]
#data_connect = data_connect[data_connect.Type .!= "NMJ", :]#FVA. Why eliminate the neuro_motor_junctions=?
###@assert unique(data_connect[!,:Type]) == ["S", "EJ", "NMJ"]


# #FVA: do it with a dictionary
# data_connect =
#     @chain data_connect begin
#         @mutate( IndexSending = fIndexByName(Sending),
#                  indexReceiving = fIndexByName(Receiving)
#             )
#     end;
# FVA: this code is idiosyncratic to VB and does not use dplyr-like processing
# # Create a new dictionary to append the indexes
# data_index = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# # Create two vectors to store the indeces before appending
# sending_value = Vector{Int64}()
# receiving_value = Vector{Int64}()
# i = 1
# # Iterate through each row and append to the new dataframe with the sending and receiving value
# for row in eachrow(data_connect)
#     for value in eachrow(neuron_position_sorted)
#         if row[:Sending] == value[:Neuron]
#             append!(sending_value, value[:Index])
#         end
#         if row[:Receiving] == value[:Neuron]
#             append!(receiving_value, value[:Index])
#         end    
#     end
#     push!(data_index, (sending_value[i], receiving_value[i]))
#     global i = i + 1
# end
# # Concatenate horizontally the dictionary of data_connect and indices
# data_connect = hcat(data_connect, data_index)

# Select the gap juntions and join them
#data_connect_gap = @filter(data_connect, Type == "EJ");
# data_connect_G = data_connect[data_connect.Type .== "G", :]
# data_connect_EJ = data_connect[data_connect.Type .== "EJ", :]
# data_connect_gap = vcat(data_connect_G, data_connect_EJ)

# # Select only the synaptic connections and join them
# data_connect_synaptic = @filter(data_connect, Type == "S");
# # data_connect_S = data_connect[data_connect.Type .== "S", :]
# # data_connect_Sp = data_connect[data_connect.Type .== "Sp", :]
# # data_connect_synaptic = vcat(data_connect_S, data_connect_Sp)
# println(describe(data_connect_gap))
# println(describe(data_connect_synaptic))
# println("Done!")

# #%% Dataframes and cleaning
# print("1.2 TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES...")
# # MONOAMINES AND NEUROPEPTIDES
# # From all the dataframe of neuron_pos_neuron only select the first three columns
# #data_connect_monoamine = data_connect_monoamine[:, 1:3]
# data_connect_monoamine =
#     @chain data_connect_monoamine begin
#         @rename(Sending = Neuron1, Receiving = Neuron2)
#         @mutate(IndexSending = fIndexByName(Sending),
#                 IndexReceiving = fIndexByName(Receiving))
#     end;
# println(describe(data_connect_monoamine))
# println("Done!")

#%% TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES
# # Create a new dictionary to append the indexes
# data_index_mono = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# # Create two vectors to store the indeces before appending
# sending_value_mono = Vector{Int64}()
# receiving_value_mono = Vector{Int64}()
# i = 1
# for row in eachrow(data_connect_monoamine)
#     for value in eachrow(neuron_position_sorted)
#         if row[:Neuron1] == value[:Neuron]
#             append!(sending_value_mono, value[:Index])
#         end
#         if row[:Neuron2] == value[:Neuron]
#             append!(receiving_value_mono, value[:Index])
#         end    
#     end
#     push!(data_index_mono, (sending_value_mono[i], receiving_value_mono[i]))
#     global i = i+1
# end
# # Concatenate horizontally the dictionary of data_connect and indeces
# data_connect_monoamine = hcat(data_connect_monoamine, data_index_mono)


# print("1.3 TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR neuropeptides...")
# data_connect_neuropep =
#     @chain data_connect_neuropep begin
#         @rename(Sending = Neuron1, Receiving = Neuron2)
#         @mutate(IndexSending = fIndexByName(Sending),
#                 IndexReceiving = fIndexByName(Receiving))
#     end;
# println(describe(data_connect_neuropep))
# println("Done!")

# # Create a new dictionary to append the indexes
# data_index_neuropep = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# # Create two vectors to store the indeces before appending
# sending_value_neuropep = Vector{Int64}()
# receiving_value_neuropep = Vector{Int64}()
# i = 1
# for row in eachrow(data_connect_neuropep)
#     for value in eachrow(neuron_position_sorted)
#         if row[:Neuron1] == value[:Neuron]
#             append!(sending_value_neuropep, value[:Index])
#         end
#         if row[:Neuron2] == value[:Neuron]
#             append!(receiving_value_neuropep, value[:Index])
#         end    
#     end
#     push!(data_index_neuropep, (sending_value_neuropep[i], receiving_value_neuropep[i]))
#     global i = i+1
# end
# # Concatenate horizontally the dictionary of data_connect and indeces
# data_connect_neuropep = hcat(data_connect_neuropep, data_index_neuropep)

# print("1.5. Table of neurotransmitters and deduced
# inhibitory/excitatory activity")
# #exploring the neurotransmitters to decide the inhibitory character of
# #neuron
# names(neuron_by_neurotransmitter)
# unique(neuron_by_neurotransmitter.Neurotransmitter1)
# unique(neuron_by_neurotransmitter.Neurotransmitter2)
# neuron_by_neurotransmitter =
#     @chain neuron_by_neurotransmitter begin
#        # @select(-(IndexNeuron))
#     #end    
#         #@mutate(IndexNeuron=fIndexByName(Neuron))
#         transform(:Neuron => ByRow(fIndexByName) => :IndexNeuron)
#         transform(
#             [:Neurotransmitter1, :Neurotransmitter2] =>
#                       ByRow((n1, n2) -> str_detect(string(n1, " ", n2), "GABA")) =>
#                       :Inhibitory,
#             [:Neurotransmitter1, :Neurotransmitter2] =>
#                       ByRow((n1, n2) ->
#                          str_detect(string(n1, " ", n2), "Glutamate|choline"))
#                   => :Excitatory
#                   )
#         #        Inhibitory=str_detect.(Neurotransmitter1, "GABA")
#     end
# println("NEURON DATA BY NEUROTRANSMITTER")
# println(describe(neuron_by_neurotransmitter))
# println("Done!")

" IMPORTANT VARIABLES:
data_connect
data_connect_monoamine
data_connect_neuropep
neuron_list
"


# Conceptually all the data should be saved to the processed directory. 
print("Saving all the preprocessed data...")
CSV.write(datadir("exp_pro", "data_connect.csv"), data_connect)
#CSV.write(datadir("exp_pro", "data_connect_gap.csv"), data_connect_gap)
CSV.write(datadir("exp_pro", "data_connect_monoamine.csv"), data_connect_monoamine)
CSV.write(datadir("exp_pro", "data_connect_neuropeptide.csv"), data_connect_neuropep)
#CSV.write(datadir("exp_pro", "neuron_position_sorted.csv"), neuron_position_sorted)
#CSV.write(datadir("exp_pro", "neuron_by_neurotransmitter.csv"), neuron_by_neurotransmitter)
CSV.write(datadir("exp_pro", "neuron_list.csv"), neuron_list)
println("Done!")


println("4. Environment description")
using Pkg;Pkg.status()


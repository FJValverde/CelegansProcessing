#! /usr/local/bin/julia
using DrWatson
@quickactivate "Celegans"

#%% Import packages
using CSV
using DataFrames
using XLSX

#%% Ignored packages for data loading and preprocessing
# using ColorSchemes
# using Colors
# using DiffEqBase
# using DifferentialEquations
# using LinearAlgebra
# using ModelingToolkit

# FVA: so dar 07/08/23

#%% DATA IMPORT
"The following lines extract the data from different files and stores as a dataframe"

data_connect_phar = CSV.read(datadir("exp_raw","ConexionsPharyngeal.csv"), DataFrame) 
# data_connect_phar = CSV.read("RawData/ConexionsPharyngeal.csv", DataFrame)
"""
Sending: Name of sending neuron
Receiving: Name of receiving neuron
Number: Number of synapses between the given neuron pair.
Type: Type of synapse: S: synaptic; G: gap.
"""


data_connect_neuron = DataFrame(XLSX.readtable(datadir("exp_raw","NeuronConnect.xlsx"), "Sheet1"))
# data_connect_neuron = DataFrame(XLSX.readtable("RawData/NeuronConnect.xlsx", "Sheet1"))
"""
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

# data_type_neuron = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet1"))
data_type_neuron = DataFrame(XLSX.readtable(datadir("exp_raw","NeuronType.xlsx"), "Sheet1"))
"""
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

                             
                                     
data_type_phar = DataFrame(XLSX.readtable(datadir("exp_raw","NeuronType.xlsx"), "Sheet2"))
# data_type_phar = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet2"))
"""
Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
"""

data_connect_monoamine = CSV.read(datadir("exp_raw","MonoaminesConnect.csv"), DataFrame)
# data_connect_monoamine = CSV.read("RawData/MonoaminesConnect.csv", DataFrame)
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type: Name of the monoamine
Monoamine: MA stands for monoamine
Specific: The specific type of monoamine
"""


data_connect_neuropep = CSV.read(datadir("exp_raw","NeuropeptidesConnect.csv"), DataFrame)
# data_connect_neuropep = CSV.read("RawData/NeuropeptidesConnect.csv", DataFrame)
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type1: One neuropeptide
Type2: Another neuropeptide
"""


neurotransmitter = DataFrame(XLSX.readtable(datadir("exp_raw","Neurotransmitters.xlsx"), "Sheet1"))
# neurotransmitter = DataFrame(XLSX.readtable("RawData/Neurotransmitters.xlsx", "Sheet1"))
"""These tables originate from a set of tables found in Pereira et al., 2015 
(http://elifesciences.org/content/4/e12432v2), combined with information compiled in Loer & Rand,
2016 & 2022, WormAtlas, and  from Gendrel et al., 2016 and Serrano-Saiz et al., 2017.						

Neuron: Name of the neuron
Neurotransmitter1: Main neurotransmitter
Neurotransmitter2: If the nueron sends another neurotransmitter it is stated here
"""


""" IMPORTANT VARIABLES: the DataFrames read in
1. data_connect_phar
2. data_connect_neuron
3. data_type_neuron
4. data_type_phar
5. data_connect_monoamine
6. data_connect_neuropep
7. neurotransmitter"""



#%% DATA FRAMES AND CLEANING
println("DATAFRAMES CLEANING")
# GAP AND SYNAPTIC
# From all the dataframe of data_type_neuron only select the first three columns
data_type_neuron = data_type_neuron[:, 1:3]
# Concatenate the two dataframes about neuron type
data_type = vcat(data_type_phar, data_type_neuron)
# Add the two neurons that do not have any connections
push!(data_type,["CANL", 0.61, "M"])
push!(data_type,["CANR", 0.61, "M"])
# Sort the dataframe by the position of the neuron in the soma
data_type_sort_soma = sort!(data_type, [:"Soma Position"])
# FVA: maybe it would be advisable to sort here also wrt distance from the tip?
# Q: Is this what the SomaPosition means roughly                            
# Add a column to know for the future the number of the neuron
data_type_sort_soma.Place = 1:302

" IMPORTANT VARIABLES:
data_type_sort_soma #FVA: As the variable gathering all GJ and Synapsis 
"

# Interchange the "Number" and "Type" columns
select!(data_connect_phar, [:Sending,:Receiving,:Type, :Number])
# Rename columns so they match the names of "data_connect_phar"
rename!(data_connect_neuron,[:Sending,:Receiving,:Type, :Number])
# Concatenate the two dataframes about neuron connectivity
data_connect = vcat(data_connect_phar, data_connect_neuron)

# Eliminate the rows in which the "Type" is either receive, receive-poly or NMJ.
data_connect = data_connect[data_connect.Type .!= "R", :]
data_connect = data_connect[data_connect.Type .!= "Rp", :]
data_connect = data_connect[data_connect.Type .!= "NMJ", :]

#%% Data homogeneization
println("TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR GAP AND SYNAPTIC CONNECTIONS.")
# Create a new dictionary to append the indexes
data_index = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# Create two vectors to store the indeces before appending
sending_value = Vector{Int64}()
receiving_value = Vector{Int64}()
i = 1
# Iterate through each row and append to the new dataframe with the sending and receiving value
for row in eachrow(data_connect)
    for value in eachrow(data_type_sort_soma)
        if row[:Sending] == value[:Neuron]
            append!(sending_value, value[:Place])
        end
        if row[:Receiving] == value[:Neuron]
            append!(receiving_value, value[:Place])
        end    
    end
    push!(data_index, (sending_value[i], receiving_value[i]))
    global i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indices
data_connect = hcat(data_connect, data_index)

# Select the gap juntions and join them
data_connect_G = data_connect[data_connect.Type .== "G", :]
data_connect_EJ = data_connect[data_connect.Type .== "EJ", :]
data_connect_gap = vcat(data_connect_G, data_connect_EJ)

# Select only the synaptic connections and join them
data_connect_S = data_connect[data_connect.Type .== "S", :]
data_connect_Sp = data_connect[data_connect.Type .== "Sp", :]
data_connect_synaptic = vcat(data_connect_S, data_connect_Sp)

#%% Dataframes and cleaning
println("FURTHER DATAFRAMES CLEANING")
# MONOAMINES AND NUEROPEPTIDES
# From all the dataframe of data_type_neuron only select the first three columns
data_connect_monoamine = data_connect_monoamine[:, 1:3]

#%% TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES
# Create a new dictionary to append the indexes
data_index_mono = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# Create two vectors to store the indeces before appending
sending_value_mono = Vector{Int64}()
receiving_value_mono = Vector{Int64}()
i = 1
for row in eachrow(data_connect_monoamine)
    for value in eachrow(data_type_sort_soma)
        if row[:Neuron1] == value[:Neuron]
            append!(sending_value_mono, value[:Place])
        end
        if row[:Neuron2] == value[:Neuron]
            append!(receiving_value_mono, value[:Place])
        end    
    end
    push!(data_index_mono, (sending_value_mono[i], receiving_value_mono[i]))
    global i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect_monoamine = hcat(data_connect_monoamine, data_index_mono)



"TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES"
# Create a new dictionary to append the indexes
data_index_neuropep = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# Create two vectors to store the indeces before appending
sending_value_neuropep = Vector{Int64}()
receiving_value_neuropep = Vector{Int64}()
i = 1
for row in eachrow(data_connect_neuropep)
    for value in eachrow(data_type_sort_soma)
        if row[:Neuron1] == value[:Neuron]
            append!(sending_value_neuropep, value[:Place])
        end
        if row[:Neuron2] == value[:Neuron]
            append!(receiving_value_neuropep, value[:Place])
        end    
    end
    push!(data_index_neuropep, (sending_value_neuropep[i], receiving_value_neuropep[i]))
    global i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect_neuropep = hcat(data_connect_neuropep, data_index_neuropep)


" IMPORTANT VARIABLES:
data_connect_synaptic
data_connect_gap
data_connect_monoamine
data_connect_neuropep
"


# Conceptually all the data should be saved to the processed directory. 
print("Saving all the preprocessed data...")
CSV.write(datadir("exp_pro", "data_connect_synaptic.csv"), data_connect_synaptic)
CSV.write(datadir("exp_pro", "data_connect_gap.csv"), data_connect_gap)
CSV.write(datadir("exp_pro", "data_connect_monoamine.csv"), data_connect_monoamine)
CSV.write(datadir("exp_pro", "data_connect_neuropeptide.csv"), data_connect_neuropep)
println("Done!")


println("4. Environment description")
using Pkg;Pkg.status()


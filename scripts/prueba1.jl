#! /usr/local/bin/julia

#%% Building the environment
using DataFrames
using CSV
using XLSX
using LinearAlgebra
using ColorSchemes
using Colors
using Plots
using DifferentialEquations
using ModelingToolkit
using DiffEqBase
using DynamicalSystems
using InteractiveChaos
using CairoMakie

#using PlotlyBase
#plotly() # FVA: the least intrusive backend, but needs the previous line

## %% Reading the data in
println("The following lines extract the data from different files and stores as a dataframe")


data_connect_phar = CSV.read("RawData/ConexionsPharyngeal.csv", DataFrame)
"""
Sending: Name of sending neuron
Receiving: Name of receiving neuron
Number: Number of synapses between the given neuron pair.
Type: Type of synapse: S: synaptic; G: gap.
"""


data_connect_neuron = DataFrame(XLSX.readtable("RawData/NeuronConnect.xlsx", "Sheet1"))
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

data_type_neuron = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet1"))
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


data_type_phar = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet2"))
"""
Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
"""

data_connect_monoamine = CSV.read("RawData/MonoaminesConnect.csv", DataFrame)
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type: Name of the monoamine
Monoamine: MA stands for monoamine
Specific: The specific type of monoamine
"""


data_connect_neuropep = CSV.read("RawData/NeuropeptidesConnect.csv", DataFrame)
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type1: One neuropeptide
Type2: Another neuropeptide
"""


neurotransmitter = DataFrame(XLSX.readtable("RawData/Neurotransmitters.xlsx", "Sheet1"))
"""These tables originate from a set of tables found in Pereira et al., 2015 
(http://elifesciences.org/content/4/e12432v2), combined with information compiled in Loer & Rand,
2016 & 2022, WormAtlas, and  from Gendrel et al., 2016 and Serrano-Saiz et al., 2017.						

Neuron: Name of the neuron
Neurotransmitter1: Main neurotransmitter
Neurotransmitter2: If the nueron sends another neurotransmitter it is stated here
"""


""" IMPORTANT VARIABLES:
data_connect_phar
data_connect_neuron
data_type_neuron
data_type_phar
data_connect_monoamine
data_connect_neuropep
neurotransmitter"""



##%% Cleaning dataframes
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
# Add a column to know for the future the number of the neuron
data_type_sort_soma.Place = 1:302

println(" IMPORTANT VARIABLES:
data_type_sort_soma
")

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


println("TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR GAP AND SYNAPTIC")
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
    i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect = hcat(data_connect, data_index)

# Select only the synaptic connections and join them
data_connect_S = data_connect[data_connect.Type .== "S", :]
data_connect_Sp = data_connect[data_connect.Type .== "Sp", :]
data_connect_synaptic = vcat(data_connect_S, data_connect_Sp)

# Select the gap juntions and join them
data_connect_G = data_connect[data_connect.Type .== "G", :]
data_connect_EJ = data_connect[data_connect.Type .== "EJ", :]
data_connect_gap = vcat(data_connect_G, data_connect_EJ)



println("DATAFRAMES CLEANING")
# MONOAMINES AND NUEROPEPTIDES
# From all the dataframe of data_type_neuron only select the first three columns
data_connect_monoamine = data_connect_monoamine[:, 1:3]

# TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES
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
    i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect_monoamine = hcat(data_connect_monoamine, data_index_mono)



println("TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES")
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
    i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect_neuropep = hcat(data_connect_neuropep, data_index_neuropep)


" IMPORTANT VARIABLES:
data_connect_synaptic
data_connect_gap
data_connect_monoamine
data_connect_neuropep
" 


println("Plot the number of connections of the synaptic and gap")
# Create two matrices for storing the connectivity information
synaptic_number = zeros(302, 302)           # To store the number of connections
synaptic_connections = zeros(302, 302)      # To store if there is a connection
# Iterate through the datatframe and add the information to the corresponding point
for link in eachrow(data_connect_synaptic)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    synaptic_number[from_index, to_index] = value
    synaptic_connections[from_index, to_index] = 1
end
# Plot both matrices
Plots.spy(synaptic_number, plot_title= "Number of synaptic connections", xlabel = "Sending neuron index", ylabel = "Rceiving neuron index", color=:berlin)
Plots.spy(synaptic_connections, plot_title= "Synaptic connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")


# For viewing the information of the frequency of the number of connections
data_s = vec(synaptic_number)           # Flatten the matrix into a 1D array
data_s = filter(x -> x != 0, data_s)    # Take out the values that are 0 for a true result
# Generate histogram
histogram(data_s, bins=:scott, 
    xticks = unique(data_s), legend = false, normalize=:probability,
    xlabel = "Number of connections between two neurons", ylabel = "Frequency", title = "Synaptic frequency histogram")
histogram!(size=(760,500))



# Repeat the previous process for GAP junctions
gap_number = zeros(302, 302)
gap_connections = zeros(302,302)
for link in eachrow(data_connect_gap)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    gap_number[from_index, to_index] = value
    gap_connections[from_index, to_index] = 1
end

Plots.spy(gap_number, plot_title= "Number of gap connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index", color=:berlin)
Plots.spy(gap_connections, plot_title= "Gap connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")

data = vec(gap_number)
data = filter(x -> x != 0, data)
# Generate histogram
histogram(data, bins=:scott, 
    xticks = unique(data), legend = false, normalize=:probability,
    xlabel = "Number of connections between two neurons", ylabel = "Frequency", title = "Gap frecuency histogram")
histogram!(size=(800,500))



"""ESTO ES EL SURPRISAL (NO ESTA BIEN TODAVIA)"""
# Calculate surprisal values for the data
surprisal_data = -log2.(data)

# Sort the data and surprisal values
sorted_data, sorted_surprisal = sort(data), surprisal_data[sortperm(data)]

# Plot the surprisal curve
Plots.plot!(sorted_data, sorted_surprisal, 
    linecolor=:blue, linewidth=2, 
    xlabel="Values", ylabel="Surprisal", title="Surprisal Curve")


#gap_number = log.(gap_number .+ 0.7)

"AQUI TERMINA EL SURPRISAL"


"Plot the connections by Monoamines"
# Create a matrix for each type
mono_connections_tyr = zeros(302,302)
mono_connections_oct = zeros(302,302)
mono_connections_dop = zeros(302,302)
mono_connections_ser = zeros(302,302)
# Add each connection to the corresponding matrix
for link in eachrow(data_connect_monoamine)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    if link[:Type] == "tyramine"
        mono_connections_tyr[from_index, to_index] = 1
    elseif link[:Type] == "octopamine"
        mono_connections_oct[from_index, to_index] = 1
    elseif link[:Type] == "dopamine"
        mono_connections_dop[from_index, to_index] = 1
    elseif link[:Type] == "serotonin"
        mono_connections_ser[from_index, to_index] = 1
    end
end
# Plot all in the same graph
Plots.spy(mono_connections_tyr, color = :orange, plot_title= "Monoamine connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")
Plots.spy!(mono_connections_oct, color = :red)
Plots.spy!(mono_connections_dop, color = :green)
Plots.spy!(mono_connections_ser, color = :blue,legend=:outertopright)




"Plot the connections by Neuropeptides (FALTA SABER QUE SIGNIFICA CADA COLUMNA)"
# Create one matrix for each column of the dataframe
neuropep_connections_Type1 = zeros(302,302)
neuropep_connections_Type2 = zeros(302,302)
# Append the value to the matrix
for link in eachrow(data_connect_neuropep)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    if link[:Type1] !== nothing
        neuropep_connections_Type1[from_index, to_index] = 1
    end
    if link[:Type2] !== nothing
        neuropep_connections_Type2[from_index, to_index] = 1
    end
end
Plots.spy(neuropep_connections_Type1, color = :blue, plot_title= "Neuropeptides connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
Plots.spy!(neuropep_connections_Type2, color = :red)


###################################################################################
###################################################################################
###################################################################################
println("Here starts the modelling, so conceptually it is no longer an exploration.")
###################################################################################


###################################################################################
println("Kunert model")
###################################################################################

# Set variables
N = 302 # number of neurons
E_rev = -65.0 #V Reversal potential of the synapse (-65 mV)
spikeThresh = 0 #V Spiking threshold
specific_capacitance = 1 #uF/cm2
intracellular_resistivity = 0.03 #kΩ*cm
g = 100 #pS Conductance gap and synaptic (Varshney et al., 2011)
Gc = 10 #pS Cell membrane conductance 
C = 1.5 #pF Membrane capacitance (Varshney et al., 2011)
Ecell = -35.0 #mV Leakage potential 
beta = 0.125  #mV−1 constant
#C = 0.015 # Cell Membrane Capacitance
ar = 1/1.5 # Synaptic activity rise time
ad = 5/1.5 # Synaptic activity decay time


Ej = fill(-48,302) #Reversal potential E_j = 0mV for excitatory synapses and −48 mV for inhibitory synapses (Wicks et al., 1996).
Iext = fill(0,301)
pushfirst!(Iext,20)

"""
Threshold potential for each neuron is computed by imposing dVi/dt=0 (Equation 2 for C. elegans) and solving for Vi. This is equivalent to Solving the following system of linear equations

where the solution x is N × 1 vector with each entry being the threshold potential Vthreshold for the ith neuron.
M1 is a matrix of size N × N where N is the number of neurons (279 for C. elegans) with its diagonal terms populated with −Gc (cell membrane capacitance).
M2 is a diagonal matrix where diagonal term in ith row corresponds to −∑jGgij i.e., the sum of total conductivity of gap junctions for the ith neuron.
M3 is a diagonal matrix where its ith diagonal term corresponds to −∑jseqGsij, where seq=arar+2ad and Gsij is maximum total conductivity of synapses to i from j. Note that seq is obtained by imposing dsidt=0 and synaptic activation Φ = 1/2 in Equation 5.
b1=Gc∗Ec where Ec is a 1D vector of size N × 1 in which all of its elements are Ec (leakage potential).
b3=Gs⋅(s∗eqEj) where Ej is a 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory).
Iext is the input stimuli vector where its ith element determines the input current amplitude for the ith neuron.
"""




println("Computing of the threshold potential (Vth)")
# Ax = b ===> A = M1+M2+M3; b = -b1-b3-Iext  (Equations 11 and 12 from Kunert 2014)

println("Computing the multilayered model by accumulating connectomes")
# M1 computation: GAP junctions
Gc_diag = zeros(302,302)    # Creation of a matrix to append the values
for i in 1:N                # In the diagonal append the constant of membrane conductance
    Gc_diag[i, i] = Gc
end
M1 = -Gc_diag               # Diagonal matrix of size N × N (N=302) with minus the membrane conductance

# M2 computation: Synaptic connections
gap_number_cond = g * gap_number       # Multiply the conductance of each channel (g=100pS) by the number of connections
gap_diag_array = sum(gap_number_cond, dims = 2)   # Sum of total conductivity for each neuron in its row
gap_diag = Diagonal(vec(gap_diag_array))  # Diagonal matrix with the values of above
M2 = -gap_diag  # Diagonal matrix of size N × N where the term corresponds to the sum of total conductivity of gap junctions for the ith neuron.


# M3 computation: Monoamines.
s_eq = ar / (ar + 2 * ad)       # Formula to compute the synaptic activity in equilibrium (is obtained by imposing dsi/dt=0 and synaptic activation Φ = 1/2 in Kunert equations)
synaptic_number_cond = g * synaptic_number  # Compute the maximum total conductivity of synapses to i from j
mult = s_eq * synaptic_number_cond          # Multiply both terms 
mult_diag_array = sum(mult, dims = 2)       # Sum of total conductivity multiplied by Seq for each neuron in its row
mult_diag = Diagonal(vec(mult_diag_array))
M3 = -mult_diag                             # Diagonal matrix of size N × N 

# b1 computation
Ec = fill(Ecell, N)     # 1D vector of size N × 1 in which all of its elements are Ec (leakage potential)
b1 = Gc * Ec            # Cell membrane conductance * Ec

# b3 computation
#Ej -> 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory)
b3 = g .* (s_eq * Ej)   # 1D vector of size N × 1 

# Compute A and b
A = M1 + M2 + M3
b = - b1 - b3 - Iext

# Solve Ax=b for x
A_inv = inv(A)      # Inverse of a
Vth = A_inv * b     # Vth is a N × 1 vector that contains the threshold potential value for each neuron
Plots.plot(Vth)           # Plot the vector Vth



function kunert_eq(u, p)
    Vi, si = u
    Gg, Gs, C, ar, ad, beta, Vth, Gc, Ecell, Ej = p
    # Parameter explanation:    Gg: matrix of N × N, corresponds to the total conductivity of gap junctions between i and j
    #                           Gs: matrix of N × N, corresponds to the maximum total conductivity of synapses between i and j
    #                           C: constant scalar, cell membrane capacitance
    #                           ar: constant scalar, correspond to the synaptic activity rise time
    #                           ad: constant scalar, correspond to the synaptic activity decay time
    #                           beta: constant scalar, width of the sigmoid function Φ (equation 6 of Kunert 2014)
    #                           Vth: vector of N × 1, corresponds to the threshold potential  
    #                           Gc: constant scalar, corresponds to the cell membrane conductance
    #                           Ecell: constant scalar, corresponds to the leakage potential
    #                           Ej: vector of N × 1, corresponds to the directionality of each neuron


    Ni = length(Vi)
    # Corresponds to equation 3 of Kunert 2014
    eq_IiGap = sum(Gg * (Vi[1] - (-70)), dims = 2)     # Sum of total conductivity multiplied by the voltage difference for each neuron in its row
    eq_IiGap = vec(eq_IiGap)                        # Vector of dimension N × 1
    
    # Corresponds to equation 4 of Kunert 2014
    eq_IiSyn = sum(Gs * si * (Vi[1] - Ej[1]), dims = 2)    # Sum of total conductivity multiplied by the voltage difference for each neuron in its row
    eq_IiSyn = vec(eq_IiSyn)                            # Vector of dimension N × 1 

    "Intento solo para la primera neurona"
    #for i in 1:Ni
    du = (-Gc * (x[1] - Ecell) - eq_IiGap[1] - eq_IiSyn[1] + 20 ) / C          # Corresponds to equation 2 of Kunert 2014
    di = ar * (1 / (1 + exp(-beta * (x[1] - Vth[1])))) * (1 - x[2]) - ad * x[2]    # Corresponds to equation 5 and 6 of Kunert 2014
    #end
    # Vi == x[1], si == x[2]
    return SVector(du, di)
end

#FVA: the sentence below is redundant
#using DifferentialEquations

# Create a vector of 302 neurons, each with two states (u and i).
neurons = zeros(302, 2)

# Define the derivatives of the states with respect to time for each neuron.
function derivatives(t, neurons)
  du = (-Gc * (neurons[:, 0] - Ecell) - eq_IiGap[1] - eq_IiSyn[1] + 20 ) / C
  di = ar * (1 / (1 + exp(-beta * (neurons[:, 0] - Vth[1])))) * (1 - neurons[:, 2]) - ad * neurons[:, 2]
  return du, di
end

#FVA: This is the line that generates the error.
# Solve the system of differential equations.
tspan = (0, 100)
solution = solve(derivatives, neurons, tspan)

# Plot the results of the simulation.
plot(tspan, solution[:, 0])
plot(tspan, solution[:, 1])

vector = [ones(10); zeros(15)]

#For the initial condition of the membrane voltages V and synaptic activity variable s, we sample the normal distribution of μ = 0 and σ = 0.94 with size 279 * 2 (for both V and s) and multiply by 10−4
# Set initial conditions (No estan bien pero por poner algo)
Vi0 = zeros(N)
si0 = fill(0.5,302)
#u0 = [Vi0; si0]
u0 = [-40, 0.5]

# Set parameters
p = [gap_number_cond, synaptic_number_cond, C, ar, ad, beta, Vth, Gc, Ecell, Ej]

# Set time span
tspan = (0.0, 0.05)  # Adjust the time span as needed
#u_begin = fill(0,604)

# Solve the differential equation
prob = ContinuousDynamicalSystem(kunert_eq, u0, p)
sol = trajectory(prob, tspan)

# Plot the solution
plot(sol, vars=(1,2))



T = 100
samp = 0.01



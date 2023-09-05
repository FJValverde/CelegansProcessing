#! /usr/local/bin/julia
#
# Author: FVA, based on V. Bokach's long monolithic script.
# 
# This script builds a Complex Dynamical System with the multilayered
#connectome of the worm C.elegans.
#
# TODO: transform this into a jupyter or Quarto document.
using DrWatson
@quickactivate "CelegansProcessing"

#%% Building the environment
#using DataFrames
#using CSV
using JLD2#load
# using DynamicalSystems
# using XLSX
using LinearAlgebra
# using ColorSchemes
# using Colors
using Plots
# using DifferentialEquations
# using ModelingToolkit
# using DiffEqBase
# using InteractiveChaos
# using CairoMakie

using Revise
#This is to be able to modify module Celegans and reload it.
# CAVEAT: use it only during Celegans development.
using Celegans

############################################################
print("1. Loading the multilayered connectome...")
connectomes = load(datadir("exp_pro",Celegans.Files.mlConnectome))
#mlConnectome = load(datadir("exp_pro",Celegans.Files.MlConnectome))
keys(connectomes)
indexByName = connectomes["indexByName"];#Overshadows the var in Celegans
newNameByIndex = connectomes["nameByIndex"];#Overshadows the var in Celegans
data_connect_synaptic_all = connectomes["dataSynapticAll"];
data_connect_gap = connectomes["dataGapJunctions"];
neuron_list = connectomes["neuronList"];
inhConnectome = connectomes["InhibitoryConnectome"];
excConnectome = connectomes["ExcitatoryConnectome"];
gapJunctionConnectome = connectomes["gapJunctionConnectome"];
println("Done!")

# Example of using the dictionary of indices and 
findall(map(x -> str_detect(x, "BAG"), newNameByIndex))

######################################################################
# Kunert Model
######################################################################

"""
Reversal potential E_j = 0mV for excitatory synapses and −45 mV for inhibitory synapses (Wicks et al., 1996).
    """
Ej = map((x -> x == "Inhibitory" ? Celegans.E_inh : Celegans.E_exc), neuron_list.NType );
@assert sum(Ej) <= sum(neuron_list.NType .== "Inhibitory")*Celegans.E_inh + sum(neuron_list.NType .== "Excitatory")*Celegans.E_exc "Warning the percentage from inhibitory to excitatory neurons does not follow the record."

"""
ALGORITHM for computing THRESHOLDS, from Neural Interactome: Interactive Simulation of a Neuronal System, Jimin Kim, William Leahy, and Eli Shlizerman, 2019, Frontiers in Computational Neuroscience

Threshold potential for each neuron is computed by imposing dVi/dt=0 (Equation 2 for C. elegans) and solving for Vi. This is equivalent to Solving the following system of linear equations

Ax = b ===> A = M1+M2+M3;

b = -b1-b3-Iext  (Equations 11 and 12 from Kunert 2014)

where the solution x is N × 1 vector with each entry being the threshold potential Vthreshold for the ith neuron.

M1 is a matrix of size N × N where N is the number of neurons (302 for a full C. elegans; 279 from the Wicks model) with its diagonal terms populated with −Gc (cell membrane capacitance), modelling current leak through the membrane.


M2 is a diagonal matrix where the diagonal term in ith row corresponds to −∑jGgij i.e., the sum of total conductivity of gap junctions for the ith neuron.

M3 is a diagonal matrix where its ith diagonal term corresponds to −∑jseqGsij, where seq=ar/ar+2ad and Gsij is maximum total conductivity of synapses to i from j. Note that seq is obtained by imposing dsi/dt=0 and synaptic activation Φ = 1/2 in Equation 5.

b1=Gc∗Ec where Ec is a 1D vector of size N × 1 in which all of its elements are Ec (leakage potential).

b3=Gs⋅(s∗eqEj) where Ej is a 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −45 mV if inhibitory).

Iext is the input stimuli vector where its ith element determines the input current amplitude for the ith neuron.
"""
print("Computing of the threshold potential (Vth)")
# Ax = b ===> A = M1+M2+M3;

# b = -b1-b3-Iext  (Equations 11 and 12 from Kunert 2014)
# M1 computation: GAP junctions
Gc_diag = zeros(nNeurons, nNeurons)    # Creation of a matrix to append the values
M1 = Diagonal(fill(-Celegans.Gc, nNeurons))
# for i in 1:nNeurons                # In the diagonal append the constant of membrane conductance
#     Gc_diag[i, i] = Gc
# end
# M1 = -Gc_diag               # Diagonal matrix of size N × N (N=302) with minus the membrane conductance

# M2 computation: diagonal with sum of total conductivity
#M2 = Diagonal(-Celegans.gᵞ * (gapJunctionConnectome * ones(nNeurons)) )
M2 = Diagonal(-Celegans.gᵞ * vec(sum(gapJunctionConnectome, dims=2)))
@assert isless(vec(sum(gapJunctionConnectome, dims=2)) - gapJunctionConnectome * ones(nNeurons), fill(eps(Float64), nNeurons))
# gap_number_cond = g * gap_number       # Multiply the conductance of each channel (g=100pS) by the number of connections
# gap_diag_array = sum(gap_number_cond, dims = 2)   # Sum of total conductivity for each neuron in its row
# gap_diag = Diagonal(vec(gap_diag_array))  # Diagonal matrix with the values of above
# M2 = -gap_diag  # Diagonal matrix of size N × N where the term corresponds to the sum of total conductivity of gap junctions for the ith neuron.

# M3 computation:
s_eq = ar / (ar + 2 * ad)       # Formula to compute the synaptic activity in eq# uilibrium (is obtained by imposing dsi/dt=0 and synaptic activation Φ = 1/2 in Kunert equations)
# synaptic_number_cond = g * synaptic_number  # Compute the maximum total conductivity of synapses to i from j
# mult = s_eq * synaptic_number_cond          # Multiply both terms 
# mult_diag_array = sum(mult, dims = 2)       # Sum of total conductivity multiplied by Seq for each neuron in its row
# mult_diag = Diagonal(vec(mult_diag_array))
# M3 = -mult_diag                             # Diagonal matrix of size N × N 
M3 = Diagonal(- gᵞ * s_eq * vec(sum(excConnectome + inhConnectome, dims=2)));

# b1 computation
b1 = Gc * fill(Ecell, nNeurons); # 1D vector of size N × 1 in which all of its elements are Ec (leakage potential)

# b3 computation
#Ej -> 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory)
b3 = gᵞ * (s_eq * Ej)   # 1D vector of size N × 1 

"""
Excitatory current for sensory experiment on CO2 sensing on BAG in C.elegans. In mV(?)
"""
Iext = zeros(nNeurons);#zero-input condition
#Iext[1] = 20
# I like the idiom above better than the one below
# Iext = fill(0,301)
# pushfirst!(Iext,20);

# Compute A and b
A = M1 + M2 + M3;
b = - b1 - b3 - Iext;

# Solve Ax=b for x
A_inv = inv(A);      # Inverse of a
"""
Vth is a N × 1 vector that contains the threshold potential value for each neuron
"""
Vth = A_inv * b
println("Done!")
        
Plots.plot(Vth)           # Plot the vector Vth
histogram(Vth,
          bins=:scott, 
          xticks = unique(Vth), legend = false, normalize=:probability,
          xlabel = "Vth of sigmoid",
          ylabel = "Frequency",
          title = "Threshold voltage histogram")
histogram!(size=(800,500))

#SO FAR: 05/09/2023 in defining the implementation as a Complex Dynamical System.


println("Computing the multilayered model by accumulating connectomes")
###############################################################################
# Model 1 using naked Dif Equations
using DifferentialEquations
##############################################################################

# Define the derivatives of the states with respect to time for each neuron.
# - u is the neuron voltage
# - i is the neuron activation state/level
# all constants imported from Celegans.Constants
Gᵞ = gᵞ * gapJunctionConnectome#The matrix of GJ conductances
Sᵞ = gᵞ * (excConnectome + inhConnectome)#Matrix of Synaptic conductances
EJ = repeat(Ej',nNeurons,1)#Same ROW is the same
@assert all(vec(EJ[rand(1:nNeurons),:] .== EJ[rand(1:nNeurons),:]))
#typeof(EJ[rand(1:nNeurons),:] .== EJ[rand(1:nNeurons),:])
function derivatives(nStates, p, t)
    #dv = (-Gc * (nStates[:, 1] - Ecell) - eq_IiGap[1] - eq_IiSyn[1] +  ) / C
    allV = repeat(nStates[:,1],1,nNeurons)#Makes the space for it. 
    dv = (-Gc * (nStates[:, 1] .- Ecell) #- eq_IiGap[1]
          -  (sum( Gᵞ .* (allV - allV'), dims=2))
          -  ( Sᵞ .* (allV - EJ)) * nStates[:, 2]) ./ C  #- eq_IiSyn[1] +  ) / C
    di = ar * (1 ./ (1 .+ exp.(-beta * (nStates[:, 1] .- Vth)))) .* (1 .- nStates[:, 2]) .- ad * nStates[:, 2]
  return hcat(dv,di)
end

"""
In place derivative function for the easier solvers.
"""
function celegansChangeState!(dS, S, p, t)
    allV = repeat(S[:,1],1,nNeurons)#Makes the space for it. 
    dS[:,1] .= (-Gc * (S[:, 1] .- Ecell) #- eq_IiGap[1]
          -  (sum( Gᵞ .* (allV - allV'), dims=2))
          -  ( Sᵞ .* (allV - EJ)) * S[:, 2]) ./ C  #- eq_IiSyn[1] +  ) / C
    dS[:,2] .= ar * (1 ./ (1 .+ exp.(-beta * (S[:, 1] .- Vth)))) .* (1 .- S[:, 2]) .- ad * S[:, 2]
  #return hcat(dv,di)
end

#FVA: This is the line that generates the error.
# Solve the system of differential equations.

# Create a vector of nNeurons neurons, each with two states (u and i).
# First in pair is neuron voltage, second is activation
"""
nState - Neuron voltage in column 1 and activation state of synapses in 2.
"""
nStates = zeros(nNeurons, 2)
dState = ODEFunction{true}(
    derivatives;
    syms=repeat(neuron_list.Neuron, 1,2)
)
dStateIIP = ODEFunction{true}(
    celegansChangeState!;
    syms=repeat(neuron_list.Neuron, 1,2)    
)
"""
For the worm we use:
- time steps of 10ms
- time spans of 10 m = 600s
For tests, timespans of 0.5m = 30s
 """
TEST = true
#TEST = false
if !TEST
    tspan = (0.0, 100.0)
else
    tspan = (0.0, 30.0)#only for half a minute.
end
fieldnames(ODEProblem)
# gusano1 = ODEProblem(dState, nStates, tspan)#Not working
# sol1 = solve(gusano1)
gusano2 = ODEProblem(dStateIIP, nStates, tspan)#Can be simulated.
timestep = 0.02# 20ms to sample at a rate of 50Hz
"""
        algo - the algorithm for solving the ODE
In DifferentialEquations.jl, some good “go-to” choices for ODEs are:

* AutoTsit5(Rosenbrock23()) handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation.
* AutoVern7(Rodas5()) handles both stiff and non-stiff equations in a way that's efficient for high accuracy.
* Tsit5() for standard non-stiff. This is the first algorithm to try in most cases.
* BS3() for fast low accuracy non-stiff.
* Vern7() for high accuracy non-stiff.
* Rodas4() or Rodas5() for small stiff equations with Julia-defined types, events, etc.
* KenCarp4() or TRBDF2() for medium-sized (100-2000 ODEs) stiff equations
* RadauIIA5() for really high accuracy stiff equations
* QNDF() for large stiff equations
    """
algo = KenCarp4()
sol2 = solve(gusano2, algo;
             alg_hints = [:stiff],#want more accuracy in solutions
             abstol = 1e-8,#This and previous param to demand more accuracy
             reltol = 1e-8,#time tolerance
             saveat = 0.02#interpolate solutions
             )
fieldnames(ODESolution)

# Plot the results of the simulation.
# TODO: tspan must be transformed to the time scale
gr()
focusNeuronPattern = "BAG"
focusNeuronIndices =
    findall(map(x -> str_detect(x, focusNeuronPattern), newNameByIndex))
focusNeurons = newNameByIndex[focusNeuronIndices]
#sol[2, 1, :], is the timeseries for the component, which is the 2nd row and 1 column.
sol2[focusNeuronIndices,1,:]
plot(sol2.t, sol2[focusNeuronIndices,1,:]')
plot(tspan, solution2[:, 1][:,2])
#TODO: plot the solutions of the problem
vector = [ones(10); zeros(15)]


###############################################################################
using ContinuousDynamicalSystems
##############################################################################
"""
The Wicks et al model in the Kunert et al formulation
"""
function kunert_eq(u, p)
    Vi, si = u#inputs,
    Gg, Gs, C, ar, ad, beta, Vth, Gc, Ecell, Ej = p#parameters
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

#For the initial condition of the membrane voltages V and synaptic activity variable s, we sample the normal distribution of μ = 0 and σ = 0.94 with size 279 * 2 (for both V and s) and multiply by 10−4
# Set initial conditions (No estan bien pero por poner algo)
Vi0 = zeros(nNeurons)
si0 = fill(0.5,nNeurons)
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



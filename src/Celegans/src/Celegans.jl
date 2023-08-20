"""
        Celegans

A module to capture idiosincrasis of the particular organism
"""
module Celegans

using DataFrames

include("celegans_constants.jl");

include("celegans_functions.jl");

include("celegans_files.jl");#Submodule Files

export nNeurons,
    # celegans_constants.jl
    E_rev,
    spikeThresh,
    Ecell,
    beta,
    specific_capacitance,
    C,
    intracellular_resistivity,
    g,
    Gc, 
    ar,
    ad,

    #celegans_functions.jl
    connectomeToMultivaluedMatrix,
    connectomeToBinaryMatrix,

    #Paths to different datasets both raw and after preprocessing
    Files
end


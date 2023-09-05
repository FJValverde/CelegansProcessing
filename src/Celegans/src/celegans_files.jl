# Set variables for Celegans
module Files

export mlConnectome, dictionaries, neuronList

"""
        Files.mlConnectome

The name of a file containing a dictionary whose components are each one of the levels of the connectome of Celegans.

@example
using JLD2
load(Celegans.Files.mlConnectome)

"""
mlConnectome = "connectome_matrices.jld2"

"""
Celegans.Files.dictionaries

The name of a file that stores a dictionary of structures compiled for the data for the whole module, including:

- indexByName: from neuron names to rank by distance from the tip

Load with:
@example
using JLD2
dictionaries = JLD2.load(datadir("exp_pro",Celegans.Files.dictionaries))
indexByName = dictionaries["indexByName"]

@example
load(datadir("exp_pro",Celegans.Files.dictionaries); indexByName)

"""
dictionaries = "dictionaries.jld2"

"""
        neuronList

The name of the file where the neuron_list dataframe is stored as a CSV file.


@example
using Celegans
neuron_list =
    CSV.read(datadir("exp_pro",Celegans.Files.neuroList), DataFrame;
             stringtype=String
"""
neuronList = "neuron_list.csv"
end

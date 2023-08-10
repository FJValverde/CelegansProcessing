" Number of neurons in the Celegans connectome."
const nNeurons = 302
# const nNeurons = 302 # number of neurons

" mV, Reversal potential of the synapse (-65 mV)"
const E_rev = -65.0 #mV, Reversal potential of the synapse (-65 mV)

const spikeThresh = 0 #V Spiking threshold
const specific_capacitance = 1 #uF/cm2
const intracellular_resistivity = 0.03 #kΩ*cm
const g = 100 #pS, Conductance gap and synaptic (Varshney et al., 2011)
const Gc = 10 #pS, Cell membrane conductance in ps
#const C = 0.015 # Cell Membrane Capacitance
const C = 1.5 #pF, Membrane capacitance (Varshney et al., 2011)
const Ecell = -35.0 #mV, Leakage potential 
const beta = 0.125  #mV−1, const
const ar = 1/1.5 # Synaptic activity rise time
const ad = 5/1.5 # Synaptic activity decay time

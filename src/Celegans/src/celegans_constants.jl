" Number of neurons in the Celegans connectome."
const nNeurons = 302
# const nNeurons = 302 # number of neurons

"""
mV, Reversal potential of the synapse (-65 mV)

But Wicks et al, 1996, reports a private communication with R.E.Davis:
- excitatory synapses: 0mV
- inhibitory synapses: -45mV

"""
const E_exc = 0#mV
const E_inh = -45#mv 
const E_rev = -65.0 #mV, Reversal potential of the synapse (-65 mV)

const spikeThresh = 0 #V Spiking threshold?????
const specific_capacitance = 1 #uF/cm2
const intracellular_resistivity = 0.03 #kΩ*cm

"""
100pS, Conductance gap and synaptic (Varshney et al., 2011)

But Wicks et al, 1996 suggests 5nS, 10x for sensibility analysis.
"""
const gᵞ = 100 #pS#FVA OK

const Gc = 10 #pS, Cell membrane conductance in ps

"""
Membrane capacitance 1(pF) (Varshney et al., 2011)
Membrane capacitance 1.5(pF) (?????)
"""
#const C = 0.015 # Cell Membrane Capacitance
const C = 1.5 #pF, Membrane capacitance (Varshney et al., 2011, say 1pF)!!

"""
Leakage potential (-35mV, Wicks et al, 1996, from David, pers.comm.)
"""
const Ecell = -35.0 #mV, Leakage potential

"""
Beta constant(mV−1) in sigmoid function.
"""
const beta = 0.125  #mV−1, const

"""
Synaptic activity rise time ( Kunert et al, 2014)
"""
const ar = 1/1.5 # FVA. Why the division? Original 1

"""
Synaptic activity decay time ( Kunert et al, 2014)
"""
const ad = 5/1.5 # FVA: Why the division? Original 5

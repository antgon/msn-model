import numpy as np
from neuron import h


def random_synapse(section, x, synapse_type,
                   stim_start=0,
                   stim_interval=1000/18,
                   stim_noise=1,
                   stim_number=1000,
                   ampa_nmda_ratio=1,
                   conn_delay=0,
                   conn_conductance=6e-4,
                   conn_threshold=0.1
                   ):
    """
    Create a random synapse in a cell section.

    Parameters
    ----------
    section : HocObject
        The cell section to receive synaptic input.
    x : numeric, range [0, 1]
        Where on the section to place the synapse.
    synapse_type : str
        Type of synapic input, 'glut' or 'gaba'.
    ampa_nmda_ratio : numeric, default=1
        The AMPA:NMDA ratio
    stim_start : numeric, default=0
        NetStim start
    stim_interval : numeric, default=1000/18
        Mean interval between two inputs in ms.
    stim_noise : numeric, default=1
        NetStim noise.
    stim_number : numeric, default=1000
        NetStim number (of what??).
    conn_delay : numeric, default=0
        NetCon delay (?).
    conn_conductance : numeric, default=6e-4
        NetConn conductance (weight).
    conn_threshold : numeric, default=0.1
        NetCon threshold.

    Returns
    -------
    synapse, netstim, netcon : HocObject
        Connection-related NEURON hoc objects 

    Notes
    -----
    Based on the Lindroos and Hellgren Kotaleski (2020) function
    random_synapse() in common_fuctions.py from the
    [striatal_SPN_lib](https://github.com/robban80/striatal_SPN_lib).
    """
    # Create synapse.
    if synapse_type == 'glut':
        synapse = h.glutamate(x, sec=section)
        synapse.ratio = ampa_nmda_ratio
    elif synapse_type == 'gaba':
        synapse = h.gaba(x, sec=section)
    else:
        raise ValueError("synapse_type must be 'glut' or 'gaba'")

    # Create stimulus (NetStim).
    netstim = h.NetStim()
    netstim.start = stim_start
    netstim.interval = stim_interval
    netstim.noise = stim_noise
    netstim.number = stim_number

    # Connect synapse (NetCon).
    netcon = h.NetCon(netstim, synapse)
    netcon.delay = conn_delay
    netcon.weight[0] = conn_conductance
    netcon.threshold = conn_threshold

    return synapse, netstim, netcon

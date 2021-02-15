"""
Build a model of a medium spiny neuron (MSN).

author: Antonio Gonzalez
"""
import numpy as np
import pandas as pd
import neuron as nrn
from neuron import h
import matplotlib.pyplot as plt

from . import paths
from .params import ModelParameters
from .funcs import random_synapse

h.load_file('stdrun.hoc')
h.load_file('import3d.hoc')
nrn.load_mechanisms(paths['mechanisms'])


def uniform_fun(distance, args, gbar):
    p0 = args[0]
    return 10**p0 * gbar


def step_fun(distance, args, gbar):
    p0, p1, p2, p3 = args
    if (distance > p2) and (distance < p3):
        value = p0
    else:
        value = p1
    return value * gbar


def get_channel_density(distance, params):
    """
    Get ion channel density (gbar, pbar) as a function of somatic
    distance.

    Parameters
    ----------
    distance : numeric
        Distance from the soma.
    params : list
        A list with four elements: [compartment, mechanism,
        density_args, gbar]; see notes below.

    Returns
    -------
    density : numeric
        The channel density value.

    Notes
    -----
    Ion channel densities are calculated for each compartment and
    mechanism as a function of somatic distance according to sigmodial,
    step, linear, or exponential functions; see Table 2 in Lindroos et
    al. (2018) and Table 1 in Lindroos and Hellgren Kotaleski (2020).
    The parameters used for these calculations are named p0, p1, p2, and
    p3 in that same Table, and that is the convention followed in this
    function. Peak conductance density is also required; the function
    get_density_params() from params.ModelParameters returns the data
    needed for these calculations in a suitable format.

    See also
    --------
    .params.ModelParameters
    """
    compartment, mech, args, gbar = params
    if compartment == 'dend':
        if mech == 'naf':
            p0, p1, p2, p3 = args
            density = gbar * 10**p0 * (1 - p1 + p1 / (
                1 + np.exp((distance-p2)/p3)))
        elif mech == 'kaf':
            p0, p1, p2, p3 = args
            density = gbar * 10**p0 * (1 + p1/(
                1 + np.exp((distance-p2)/p3)))
        elif mech == 'kas':
            p0, p2, p3 = args
            density = gbar * 10**p0 * (0.1 + 0.9/(
                1 + np.exp((distance-p2)/p3)))
        elif (mech == 'kir') or (mech == 'sk'):
            p0 = args[0]
            density = gbar * 10**p0
        elif mech == 'can':
            p0, p1, p2, p3 = args
            density = 10**p0 * (1 - p1 + p1/(
                1 + np.exp((distance-p2)/p3)))
        elif (mech == 'cav32') or (mech == 'cav33'):
            p0, p2, p3 = args
            density = 10**p0 / (1 + np.exp((distance-p2)/p3))
        else:
            density = uniform_fun(distance, args, gbar)
    elif compartment == 'axon':
        if (mech == 'kas') or (mech == 'Im'):
            density = uniform_fun(distance, args, gbar)
        elif mech == 'naf':
            density = step_fun(distance, args, gbar)
    elif compartment == 'soma':
        density = uniform_fun(distance, args, gbar)

    # In the function `calculate_distribution` (in Lindroos's
    # MSN_builder.py) the resulting density may be negative, in which
    # case they set the density to 0 instead. I do not know why or when
    # the density may be negative. This assert is to look into this if
    # that ever happened.
    assert (density >= 0), "Channel distribution is negative"

    return density


class MSN:
    """
    Build a model of a MSN.

    Attributes
    ----------
    type : str
        Cell type, 'imsn' or 'dmsn'
    index : int
        Cell index
    distrib_params : list
        Channel distribution parameters (compartment, mechanism, args,
        gbar)
    rheobase : numeric
    v_init : numeric
        Initialisation membrane voltage

    Methods
    -------
    add_bg_noise(gaba_freq=4, glut_freq=12, syn_fact=[],
                 gaba_mod=0, dend_only=False, delays=[]
        Add background noise
    remove_bg_noise()
        Remove background noise

    Notes
    -----
    This class is based on 'MSN_builder.py', distributed with the NEURON
    model by Lindroos and Hellgren Kotaleski (2020), available from
    ModelDB (accession number 266775).
    """
    def __init__(self, cell_type, cell_index, v_init=-80):
        """
        Parameters
        ----------
        cell_type : str
            Cell type to model, one of 'dmsn' or 'imsn'.
        cell_index : int
            Cell to model, from the set of iMSN and dMSN provided by
            Lindroos et al. In that set there are n=71 dMSNs and n=34
            iMSNs.
        v_init : numeric, default=-80
            Initialisation membrane voltage.
        """
        # self._gid = gid
        self.type = cell_type
        self.index = cell_index

        # Load parameters
        params = ModelParameters()
        self.density_params = params.get_density_params(cell_type, cell_index)
        self.rheobase = params.get_rheobase(cell_type, cell_index)
        self._morphology_file = params.get_morphology_path(cell_type)
        gbar_pas = params.get_gbar(cell_type)
        self._gbar_pas = gbar_pas.value[gbar_pas.mechanism == 'pas'].values[0]

        # Create cell
        self._setup_morphology()
        self._setup_mechanisms()
        self._setup_biophysics()
        self._setup_density()
        h.celsius = 35
        self.v_init = v_init

        # Additional containers
        self._bg_noise = []

    def _setup_morphology(self):
        # Import morphology from SWC file.
        cell = h.Import3d_SWC_read()
        cell.input(self._morphology_file)
        i3d = h.Import3d_GUI(cell, 0)
        i3d.instantiate(self)
        # h.define_shape() # Is this needed? This line is in the original
        # Lindroos MSN but I am not sure what it does.

        # There should only be one soma
        assert(len(self.soma) == 1)
        self.soma = self.soma[0]

        # Spatial discretisation.
        for sec in self.all:
            if 'axon' in sec.name():
                sec.nseg = 2
            else:
                sec.nseg = 2 * int(sec.L/40) + 1

    def _setup_mechanisms(self):
        dendritic_channels = (['naf', 'kaf', 'kas', 'kir', 'sk', 'can',
                               'cav32', 'cav33', 'kdr', 'cal12', 'cal13',
                               'car', 'bk'] +
                              ['cadyn', 'caldyn'])
        somatic_channels = (['naf', 'kaf', 'kas', 'kdr', 'bk', 'cal12',
                             'cal13', 'car', 'can', 'sk', 'kir'] +
                            ['cadyn', 'caldyn'])
        axonal_channels = ['naf', 'kas', 'Im']

        for sec in self.all:
            if 'soma' in sec.name():
                mechanisms = somatic_channels
            elif 'axon' in sec.name():
                mechanisms = axonal_channels
            elif 'dend' in sec.name():
                mechanisms = dendritic_channels
            else:
                raise ValueError('Unknown section: ' + sec.name())
            for mechanism in mechanisms:
                sec.insert(mechanism)

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 150
            sec.cm = 1
            sec.insert('pas')
            for seg in sec:
                seg.pas.g = self._gbar_pas
                seg.pas.e = -70
            sec.ena = 50
            sec.ek = -85

    def _setup_density(self):
        """
        Setup ion channel density.
        """
        for params in self.density_params:
            compartment, mechanism, args, gbar = params
            if mechanism.startswith('ca'):
                prefix = 'pbar_'
            else:
                prefix = 'gbar_'
            h.distance(sec=self.soma)
            for section in self.all:
                if compartment in section.name():
                    for segment in section:
                        distance = h.distance(segment.x, sec=section)
                        density = get_channel_density(distance, params)
                        cmd = f'segment.{prefix}{mechanism} = {density}'
                        exec(cmd)

    def add_bg_noise(self, gaba_freq=4, glut_freq=12, dend_only=False,
                     ampa_scale_factor=None, nmda_scale_factor=None,
                     gaba_scale_factor=None, delays=[]):
        """
        Add background noise.

        Creates random GABA and glutamate synaptic inputs and connects these
        to the cell sections to create background noise.

        Parameters
        ----------
        gaba_freq : numeric, default=4
            Frequency of GABAergic inputs (number of inputs per s)
        glut_freq : numeric, default=12
            Frequency of glutamatergic inputs (number of input per s)
        dend_only : bool, default=False
            Whether glut and GABA inputs are applied only to dendrites or to
            all cell sections.
        ampa_scale_factor,
        nmda_scale_factor : None or numeric, default=None
            If not None, scale glutamatergic input by these factors.
            (These two were parameter `syn_fact` in the original
            Lindroos function.)
        gaba_scale_factor : None or numeric, default=None
            If not None, scale GABA input by these factors.
            (This was parameter `gabaMod` in the original Lindroos
            function.)
        delays : list, default=[]
            TODO: What are `delays`?
            Looks like it should be a list with the same number of
            elements as there are sections in the cell. Then, each
            element will set the value of `NetStim.start` for each
            section. (It is confusing that in the random_synapse()
            function there is a NetCon.delay parameter but this delays
            here seem to be unrelated.)

        Notes
        -----
        Based on the function set_bg_noise() in common_functions.py,
        from Lindroos et al. and available from
        [GitHub](https://github.com/robban80/striatal_SPN_lib).
        """
        if dend_only is True:
            sections = [sec for sec in self.all if 'dend' in sec.name()]
        else:
            sections = self.all

        # Make sure to remove previous noise first.
        self.remove_bg_noise()

        # What is gbase? Why these default values? Where do they come
        # from?
        if self.type == 'dmsn':
            gbase = 0.3e-3
        elif self.type == 'imsn':
            gbase = 0.2e-3

        for indx, sec in enumerate(sections):
            if len(delays) == 0:
                delay = 0
            else:
                delay = delays[indx]

            # Glut synapse
            synapse, netstim, netcon = random_synapse(
                sec, x=0.5, synapse_type='glut',
                stim_interval=1000/glut_freq,
                conn_conductance=gbase, stim_start=delay)
            synapse.ratio = 1
            if ampa_scale_factor:
                synapse.ampa_scale_factor = ampa_scale_factor
            if nmda_scale_factor:
                synapse.nmda_scale_factor = nmda_scale_factor
            self._bg_noise.append([synapse, netstim, netcon])

            # GABA synapse
            if gaba_scale_factor:
                conductance = gbase * 3 * gaba_scale_factor  # Why times 3?
            else:
                conductance = gbase * 5  # Why times 5?
            synapse, netstim, netcon = random_synapse(
                sec, x=0.1, synapse_type='gaba',
                stim_interval=1000/gaba_freq,
                conn_conductance=conductance,
                stim_start=delay)
            self._bg_noise.append([synapse, netstim, netcon])

    def remove_bg_noise(self):
        """
        Removes background noise from the cell.

        Notes
        -----
        The way to remove synaptic inputs from a cell model seems to be
        to just set NetCon weight to 0, as mentioned in this post:

        https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3576
        """
        if len(self._bg_noise) > 0:
            for line in self._bg_noise:
                netcon = line[2]
                netcon.weight[0] = 0
            self._bg_noise = []

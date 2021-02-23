TITLE Potassium current, M-type (I_K(M))

COMMENT
A slowly activating, non-inactivating potassium current, muscarine-
sensitive. Modelled after Adams et al (1982).

A similar implentation of this current is that by Doron et al (2017),
available on ModelDB (#231427). Note, however, that their equations used
to calculate `alpha` and `beta` differ slightly from Equations 5a and 5b
in Adams et al.

Temperature compensation in this mechanism follows that in Doron's model.
It also includes a `modulation()` function, from Lindroos and Hellgren
Kotaleski (2020) (ModelDB #266775), use to simulate dopamine and
acetylcholine modulation of this mechanism.

References
----------
Adams PR, Brown DA & Constanti A (1982). M-currents and other potassium
currents in bullfrog sympathetic neurones. J Physiol 330: 537-572.

Lindroos R & Hellgren Kotaleski J (2020). Predicting complex spikes in
striatal projection neurons of the direct pathway following
neuromodulation by acetylcholine and dopamine. Eur J Neurosci
(https://doi.org/10.1111/ejn.14891).


(2020) Antonio Gonzalez
ENDCOMMENT

NEURON {
    SUFFIX Im
    USEION k READ ek WRITE ik
    RANGE gbar, g, i
    RANGE damod, maxMod, level, max2, lev2
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S)  = (siemens)
    : F = (faraday) (kilocoulombs)
    : R = (k-mole) (joule/degC)
}

PARAMETER {
    gbar = 0.00001 (S/cm2)
    damod = 0
    maxMod = 1
    max2 = 1
    level = 0
    lev2 = 0
    q10 = 2.3  : Temperature sensitivity; from Doron's model.
    celsius (degC)
    temp = 22 (degC)  : Original temperature, from Adams et al (1982).
}

ASSIGNED {
    v (mV)
    ek (mV) 
    ik (mA/cm2)
    i (mA/cm2)
    g (S/cm2)
    tadj (1)
}

STATE {
    m (1)
}

INITIAL {    
    tadj = q10^((celsius - temp)/(10(degC)))  : Temperature-adjusting
                                              : factor. From Doron's model.
    m = minf(v)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gbar * m * modulation()
    i = g * (v - ek)
    ik = i
}

DERIVATIVE states {
    m' = (minf(v) - m) / taum(v)
}

FUNCTION minf (Vm (mV)) (1) {
    : After Equation 2 in Adams et al (1982).
    : That equation is:
    : minf = 1/(1 + exp(z * (v0 - Vm) * F/R/T))
    : where 
    :   z = 2.5 (1)
    :   v0 = -35 (mV)
    :   F = 96.480 (kilocoul/mole)
    :   R = 8.314 (joule/degC)
    :   T = 273.15 + 22(celsius)
    : In this implementation, I have simplified that equation by
    : reducing all constants to one number, km. Thus, `km` (below) is
    : 1/(z*F/R/T).
    LOCAL vm, km
    vm = -35 (mV)
    km = 10.17 (mV)
    minf = 1 / (1 + exp((vm - Vm) / km))
}

FUNCTION taum (Vm (mV)) (ms) {
    : After Equations 4, 5a and 5b in Adams et al (1982).
    : As with `minf`, I reduced the constants in the original equations
    : to the parameters `ktau1` and `ktau2`.
    LOCAL vtau, ktau1, ktau2, a, b
    vtau = -35 (mV)
    ktau1 = 20.3 (mV)
    ktau2 = -20.3 (mV)
    a = 3.3(/s) * exp((Vm - vtau)/ktau1)
    b = 3.3(/s) * exp((Vm - vtau)/ktau2)
    taum = 1/(a + b) * (1000)
    taum = taum/tadj
}

FUNCTION modulation() {
    modulation = 1 + damod * (
        (maxMod - 1) * level +
        (max2 - 1) * lev2)
    if (modulation < 0) {
        modulation = 0
    }    
}
TITLE NMDA synapse

COMMENT
This mechanism is based on Examples 10.3 to 10.5 in The NEURON Book
with additions and parameters from Lindross's `gluatamte.mod`
(ModelDB #266775). These additions include Mg++ block and the
`modulation()` function for modelling DA and Ach modulation.

This mechanism contributes (with AMPA) to ical.

NOTE: The original glutamate.mod file fails modlunit.

(2020) Antonio Gonzalez
ENDCOMMENT

NEURON {
    POINT_PROCESS nmda
    USEION cal WRITE ical VALENCE 2
    RANGE e, g, i
    RANGE damod, maxMod, level, max2, lev2
    NONSPECIFIC_CURRENT i
}

UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
    (mM) = (milli/liter)
}

PARAMETER {
    e = 0 (mV)
    q = 2
    tau1 = 5.52 (ms) : Rise time constant.
    tau2 = 231 (ms) : Decay time constant. Must be greater than tau1.
    
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
    scale_factor = 1 : Scales the total current.

    ca_ratio = 0.1 : Fraction of current across NMDA that is carried by Ca++.
}

ASSIGNED {
    v (mV) 
    i (nA)
    itotal (nA)
    g (uS)
    ical (nA)
    factor
}

STATE {
    a (uS)
    b (uS)
}

INITIAL {
    LOCAL tp
	a = 0
	b = 0
	tp = (tau1 * tau2)/(tau2 - tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	g = (b - a) * modulation()
	itotal = g * (v - e) * mgblock(v) * scale_factor
    ical = itotal * ca_ratio
    i = itotal * (1 - ca_ratio)
}

DERIVATIVE states {
	a' = -a/tau1 * q
	b' = -b/tau2 * q
}

NET_RECEIVE(weight (uS)) {
	a = a + weight * factor
	b = b + weight * factor
}

FUNCTION mgblock(Vm (mV)) (1) {
    LOCAL mg, alpha, beta
    mg = 1 (mM)
    alpha = 0.062 (/mV)
    beta = 3.57 (mM)
    mgblock = 1 / (1 + mg * exp(-alpha * Vm) / beta)
}

: In the original function (in glutamate.mod) the modulation
: variables have different names. They were changed here to match the 
: names found it the other mechanims in the model.
: m1 = maxMod
: l1 = level
: m2 = max2
: l2 = lev2
FUNCTION modulation() (1) {
    modulation = 1 + damod * (
        (maxMod - 1) * level +
        (max2 - 1) * lev2)
    if (modulation < 0) {
        modulation = 0
    }    
}
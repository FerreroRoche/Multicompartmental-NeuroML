TITLE Mod file for component: Component(id=vClamp0 type=voltageClampTriple)

COMMENT

    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.5.5
         org.neuroml.model   v1.5.5
         jLEMS               v0.9.9.2

ENDCOMMENT

NEURON {
    POINT_PROCESS vClamp0
    ELECTRODE_CURRENT i
    RANGE weight                            : property
    RANGE active                            : parameter
    RANGE delay                             : parameter
    RANGE duration                          : parameter
    RANGE conditioningVoltage               : parameter
    RANGE testingVoltage                    : parameter
    RANGE returnVoltage                     : parameter
    RANGE simpleSeriesResistance            : parameter
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    weight = 1
    active = 0 
    delay = 50 (ms)
    duration = 200 (ms)
    conditioningVoltage = -55 (mV)
    testingVoltage = -55 (mV)
    returnVoltage = -55 (mV)
    simpleSeriesResistance = 1 (Mohm)
}

STATE {
    i (nA) 
    
}

INITIAL {
    rates()
    rates() ? To ensure correct initialisation.
    
}

BREAKPOINT {
    
    rates()
    if (active  == 1 && t <  delay) {
        i = weight  * (  conditioningVoltage   - v) /  simpleSeriesResistance ? standard OnCondition
    }
    
    if (active  == 1 && t >=  delay) {
        i = weight  * (  testingVoltage   - v) /  simpleSeriesResistance ? standard OnCondition
    }
    
    if (active  == 1 && t >  duration  +  delay) {
        i = weight  * (  returnVoltage   - v) /  simpleSeriesResistance ? standard OnCondition
    }
    
    
}

PROCEDURE rates() {
    
    
     
    
}

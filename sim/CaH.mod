TITLE High Threshold Calcium
: (Jung Lee 2013)

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX CaH
        USEION ca READ eca WRITE ica
	RANGE gbar
}

PARAMETER {
	gbar  = 0.00012 (S/cm2)	:from Pospischil et al. - 2008 - Minimal Hodgkin-Huxley type models for different classes of cortical and thalamic neurons
}

ASSIGNED { 
	v (mV)
	eca (mV)
        ica (mA/cm2)
        malpha (/ms)
        mbeta (/ms)
	:ninf
	:ntau (ms)
        
}

STATE {
        m
	:n
}

BREAKPOINT {
	SOLVE states METHOD cnexp 
	ica  = gbar*m*m*(v-eca)
	:ica = gbar*n*n*n*n*(v - eca)
}

INITIAL {
        settables(v)
        m = malpha/(malpha+mbeta)
	:n = ninf
}

DERIVATIVE states {  
	settables(v)      
	m' = ((malpha*(1-m)) - (mbeta*m))
	:n' = (ninf-n)/ntau
}

UNITSOFF

PROCEDURE settables(v (mV)) {
        TABLE malpha, mbeta  FROM -100.5 TO 95.5 WITH 200   :Offset to avoid singularity at 30mV
	:TABLE ninf, ntau  FROM -100 TO 100 WITH 200

:	malpha = (1.6) / (1+exp(-0.072*(v-5)))
:	mbeta = 0.02*(v+8.9) / (exp((v+8.9)/5) - 1)

	malpha = 1.6 / (1+exp(-0.072*(v-5)))
	mbeta = 0.02*(v+8.9) / (exp((v+8.9)/5) - 1)

	:ninf = 1.0/(1.0+exp((-v-27.0)/11.5))
	:ntau = (0.25 + 4.35*exp(-fabs(v+10)/10))*1
}

UNITSON




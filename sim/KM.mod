TITLE Potassium M Current for all cells

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX KM
        USEION k READ ek WRITE ik
	RANGE gkmbar
}

PARAMETER {
	gkmbar  = 0.004 (S/cm2)
}

ASSIGNED { 
	v (mV)
	ek (mV)
        ik (mA/cm2)
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
	ik  = gkmbar*m*(v-ek)
	:ik = gkmbar*n*n*n*n*(v - ek)
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
        LOCAL  Qs
        TABLE malpha, mbeta  FROM -100.5 TO 95.5 WITH 200   :Offset to avoid singularity at 30mV
	:TABLE ninf, ntau  FROM -100 TO 100 WITH 200

	Qs = 3.209
	malpha = (0.0001*Qs*(v+30)) / (1-exp(-(v+30)/9))
	mbeta = (-0.0001*Qs*(v+30)) / (1-exp((v+30)/9))

	:ninf = 1.0/(1.0+exp((-v-27.0)/11.5))
	:ntau = (0.25 + 4.35*exp(-fabs(v+10)/10))*1
}

UNITSON

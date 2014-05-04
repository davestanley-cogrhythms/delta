TITLE hh_ex.mod NaF and KDR channels for inhibitory cells
 
COMMENT
	Modified defualt neuron hh object for inhibitory NaF and KDR channels. Taken from Lee 2013.
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hh_ex
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        RANGE gnabar, gkbar, gna, gk, m0
        GLOBAL hinf, ninf, htau, ntau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
}
 
STATE {
        h n
}
 
ASSIGNED {
        v (mV)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        m0 hinf ninf
	htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp	: Not sure whether to leave this in or not!
        gna = gnabar*m0*m0*m0*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)
}
 
 
INITIAL {
	rates(v)
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        TABLE m0, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200

UNITSOFF
	m0 = 1.0/(  1+exp(-(v+34.5)/10)  )

	hinf = 1.0/(  1+exp((v+59.4)/10.7)  )
	htau = (0.15 + 1.15/(  1.0+exp((v+33.5)/15.0)  ))*1

	ninf = 1.0/(1.0+exp((-v-29.5)/10.0))
	ntau = (0.25 + 4.35*exp(-fabs(v+10)/10))*1

}
 
UNITSON



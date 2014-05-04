
function matlab_LTS
    dt = 0.01; % ms
    endT=10000;
    Y0 = get_Y0(-60); Y0=Y0(:);
    [T YT] = ode45(@diff_LTS,(0:endT-1)*dt,Y0);
    
    figure; plot(T,YT(:,1));
end


function Ypr = diff_LTS(t,Y)
    ek = -90;
    ena = 50;
    e_pas = -67;
    C = 1.0;
    gnabar = 100; %uS/cm2
    gkbar = 80;
    g_pas = 0.1;
    i = -10;

    v=Y(1);
    h=Y(2);
    n=Y(3);
    
    m0 = 1.0/(  1+exp(-(v+38)/10)  );

	hinf = 1.0/(  1+exp((v+58.3)/6.7)  );
	htau = 0.225 + 1.125/(  1.0+exp((v+37.0)/15.0)  );

	ninf = 1.0/(1.0+exp((-v-27.0)/11.5));
	ntau = 0.25 + 4.35*exp(-abs(v+10)/10);

    gna = gnabar*m0*m0*m0*h;
	ina = gna*(v - ena);
    gk = gkbar*n*n*n*n;
	ik = gk*(v - ek);
    ipas = g_pas*(v-e_pas);
    
    vpr = -1/C*(ina + ik + ipas + i)/1;   % mV/second. Note that should divide this by 1000 to get mV/ms.
                                             % Need to then increase g or decrease C by a factor of 1000 to compensate.
    hpr = (hinf-h)/htau;            %? - I guess this is ms
    npr = (ninf-n)/ntau;            %? - I guess this is ms

    Ypr(1) = vpr;
    Ypr(2) = hpr;
    Ypr(3) = npr;
    Ypr=Ypr(:);

end


function Y0 = get_Y0 (v)

	hinf = 1.0/(  1+exp((v+58.3)/6.7)  );
	htau = 0.225 + 1.125/(  1.0+exp((v+37.0)/15.0)  );

	ninf = 1.0/(1.0+exp((-v-27.0)/11.5));
	ntau = 0.25 + 4.35*exp(-abs(v+10)/10);
    
    Y0(1) = v;
    Y0(2) = hinf;
    Y0(3) = ninf;
    
end



% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case2
mpc.version = '2';
mpc.baseMVA =  100.00;


%% bus data
%	bus_i	type	Pd	Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV	    zone	Vmax	    Vmin
mpc.bus = [
	1	    2	 400.0	 0.0	 0.0	 1.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	2	    2	 600.0	 0.0	 0.0	 2.0	 2	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg	    Qg		Qmax	  Qmin	   Vg	    mBase	status	 Pmax	  Pmin	  H	 min_su_time	min_sd_time	
mpc.gen = [
	1	 0.0	0.0	 	30.0	 -30.0	   1.00000  100.0	 1	 1100.0	  50.0	    3;
	1	 0.0	0.0	 	30.0	 -30.0	   1.00000  100.0	 1	 1100.0	  50.0	    8;
	2	 0.0	0.0	 	30.0	 -30.0	   1.00000  100.0	 1	 1100.0	  50.0	    3;
	2	 0.0	0.0	 	30.0	 -30.0	   1.00000  100.0	 1	 1100.0	  50.0	    10;

];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 200.0	 0.0	 3	   0.000000	  5.000000	   0.000000	   0.000000;
	2	 100.0	 0.0	 3	   0.000000	  15.0100	   0.000000	   0.000000;
	2	 100.0	 0.0	 3	   0.000000	  7.500000	   0.000000	   0.000000;
	2	 500.0	 0.0	 3	   0.000000	  20.000000	   0.000000	   0.000000;

];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.000281	 0.00281	 0.00712	 4000.0	 4000.0	 4000.0	 0.0	  0.0	 1	 -180.0	 180.0;
];


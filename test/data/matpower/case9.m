% used in tests of,
% - sparce SDP implementation, possible cholesky PosDefException

function mpc = case9
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs			area	Vm			Va		baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 100.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	2	 2	 100.0	 0.0	 0.0	 0.0	 2	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	3	 2	 100.0	 0.0	 0.0	 0.0	 3	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	4	 1	 100.0	 0.0	 0.0	 0.0	 1	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	5	 1	 100.0	 20.0	 0.0	 0.0	 2	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	6	 1	 100.0	 0.0	 0.0	 0.0	 3	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	7	 1	 100.0	 40.0	 0.0	 0.0	 1	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	8	 1	 100.0	 0.0	 0.0	 0.0	 2	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
	9	 1	 100.0	 40.0	 0.0	 0.0	 3	    1.00000	    0.00000	 350.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg		Qg		 Qmax	  Qmin	 Vg		 mBase	status	Pmax  Pmin	H
mpc.gen = [
	1	 0.0	 0.0	 900.0	 -900.0	 1.0	 100.0	 1	 1200.0	 100.0	5;
	2	 0.0	 0.0	 900.0	 -900.0	 1.0	 100.0	 1	 1000.0	 100.0	10;
	3	 0.0	 0.0	 900.0	 -900.0	 1.0	 100.0	 1	 800.0	 100.0	7;
];

%% generator cost data
mpc.gencost = [
	2	 0.0	 0.0	 2	   4.1	 0.0;
	2	 0.0	 0.0	 2	   7.3	 0.0;
	2	 0.0	 0.0	 2	   1.1	 0.0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 4	 0.000	 0.058	 0.000	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 5	 0.017	 0.092	 0.158	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	5	 6	 0.039	 0.170	 0.358	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 6	 0.000	 0.059	 0.000	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 7	 0.012	 0.101	 0.209	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	7	 8	 0.009	 0.072	 0.149	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	8	 2	 0.000	 0.063	 0.000	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	8	 9	 0.032	 0.161	 0.306	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
	9	 4	 0.010	 0.085	 0.176	 900.0	 900.0	 900.0	 0.0	 0.0	 1	 -30.0	 30.0;
];

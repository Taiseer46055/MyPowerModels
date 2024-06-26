% from pglib-opf - https://github.com/power-grid-lib/pglib-opf

function mpc = case30
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	   -0.00000	 132.0	 1	    1.06000	    0.94000;
	2	 2	 21.7	 12.7	 0.0	 0.0	 1	    1.03591	   -4.11149	 132.0	 1	    1.06000	    0.94000;
	3	 1	 2.4	 1.2	 0.0	 0.0	 1	    1.01502	   -6.85372	 132.0	 1	    1.06000	    0.94000;
	4	 1	 7.6	 1.6	 0.0	 0.0	 1	    1.00446	   -8.44574	 132.0	 1	    1.06000	    0.94000;
	5	 2	 94.2	 19.0	 0.0	 0.0	 1	    0.99748	  -13.13219	 132.0	 1	    1.06000	    0.94000;
	6	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00170	  -10.15671	 132.0	 1	    1.06000	    0.94000;
	7	 1	 22.8	 10.9	 0.0	 0.0	 1	    0.99208	  -11.91706	 132.0	 1	    1.06000	    0.94000;
	8	 2	 30.0	 30.0	 0.0	 0.0	 1	    1.00241	  -10.93461	 132.0	 1	    1.06000	    0.94000;
	9	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.03671	  -13.26615	 1.0	 1	    1.06000	    0.94000;
	10	 1	 5.8	 2.0	 0.0	 19.0	 1	    1.03220	  -14.89730	 33.0	 1	    1.06000	    0.94000;
	11	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	  -13.26615	 11.0	 1	    1.06000	    0.94000;
	12	 1	 11.2	 7.5	 0.0	 0.0	 1	    1.04625	  -14.18878	 33.0	 1	    1.06000	    0.94000;
	13	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	  -14.18878	 11.0	 1	    1.06000	    0.94000;
	14	 1	 6.2	 1.6	 0.0	 0.0	 1	    1.03103	  -15.09457	 33.0	 1	    1.06000	    0.94000;
	15	 1	 8.2	 2.5	 0.0	 0.0	 1	    1.02621	  -15.17834	 33.0	 1	    1.06000	    0.94000;
	16	 1	 3.5	 1.8	 0.0	 0.0	 1	    1.03259	  -14.75320	 33.0	 1	    1.06000	    0.94000;
	17	 1	 9.0	 5.8	 0.0	 0.0	 1	    1.02726	  -15.07475	 33.0	 1	    1.06000	    0.94000;
	18	 1	 3.2	 0.9	 0.0	 0.0	 1	    1.01602	  -15.79032	 33.0	 1	    1.06000	    0.94000;
	19	 1	 9.5	 3.4	 0.0	 0.0	 1	    1.01317	  -15.95831	 33.0	 1	    1.06000	    0.94000;
	20	 1	 2.2	 0.7	 0.0	 0.0	 1	    1.01714	  -15.75153	 33.0	 1	    1.06000	    0.94000;
	21	 1	 17.5	 11.2	 0.0	 0.0	 1	    1.01982	  -15.35269	 33.0	 1	    1.06000	    0.94000;
	22	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.02041	  -15.33860	 33.0	 1	    1.06000	    0.94000;
	23	 1	 3.2	 1.6	 0.0	 0.0	 1	    1.01532	  -15.56400	 33.0	 1	    1.06000	    0.94000;
	24	 1	 8.7	 6.7	 0.0	 4.3	 1	    1.00930	  -15.72597	 33.0	 1	    1.06000	    0.94000;
	25	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.00621	  -15.29138	 33.0	 1	    1.06000	    0.94000;
	26	 1	 3.5	 2.3	 0.0	 0.0	 1	    0.98833	  -15.72053	 33.0	 1	    1.06000	    0.94000;
	27	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.01294	  -14.75624	 33.0	 1	    1.06000	    0.94000;
	28	 1	 0.0	 0.0	 0.0	 0.0	 1	    0.99824	  -10.79567	 132.0	 1	    1.06000	    0.94000;
	29	 1	 2.4	 0.9	 0.0	 0.0	 1	    0.99288	  -16.01182	 33.0	 1	    1.06000	    0.94000;
	30	 1	 10.6	 1.9	 0.0	 0.0	 1	    0.98128	  -16.91371	 33.0	 1	    1.06000	    0.94000;
];

%% generator data
%	bus	Pg			Qg		 Qmax	 Qmin	 Vg			 mBase	status	Pmax	Pmin	H
mpc.gen = [
	1	 218.839	9.372	 10.0	 0.0	 1.06	 	 100.0	 1	 	784	 	0.0		12; % COW
	2	 80.05	 	24.589	 50.0	 -40.0	 1.03591	 100.0	 1	 	100	 	0.0		7; % NG
	5	 0.0	 	32.487	 40.0	 -40.0	 0.99748	 100.0	 1	 	0	 	0.0		0; % SYNC
	8	 0.0	 	40.0	 40.0	 -10.0	 1.00241	 100.0	 1		0	 	0.0		0; % SYNC
	11	 0.0	 	11.87	 24.0	 -6.0	 1.06	 	 100.0	 1	 	0	 	0.0		0; % SYNC
	13	 0.0	 	10.414	 24.0	 -6.0	 1.06	 	 100.0	 1	 	0	 	0.0		0; % SYNC
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   0.000000	   0.521378	   0.000000; % COW
	2	 0.0	 0.0	 3	   0.000000	   1.135166	   0.000000; % NG
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000; % SYNC
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000; % SYNC
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000; % SYNC
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000; % SYNC
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.0192	 0.0575	 0.0528	 138	 138	 138	 0.0	 0.0	 1	 -30.0	 30.0;
	1	 3	 0.0452	 0.1652	 0.0408	 152	 152	 152	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 4	 0.057	 0.1737	 0.0368	 139	 139	 139	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 4	 0.0132	 0.0379	 0.0084	 135	 135	 135	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 5	 0.0472	 0.1983	 0.0418	 144	 144	 144	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 6	 0.0581	 0.1763	 0.0374	 139	 139	 139	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 6	 0.0119	 0.0414	 0.009	 148	 148	 148	 0.0	 0.0	 1	 -30.0	 30.0;
	5	 7	 0.046	 0.116	 0.0204	 127	 127	 127	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 7	 0.0267	 0.082	 0.017	 140	 140	 140	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 8	 0.012	 0.042	 0.009	 148	 148	 148	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 9	 0.0	 0.208	 0.0	 142	 142	 142	 0.978	 0.0	 1	 -30.0	 30.0;
	6	 10	 0.0	 0.556	 0.0	 53	 53	 53	 0.969	 0.0	 1	 -30.0	 30.0;
	9	 11	 0.0	 0.208	 0.0	 142	 142	 142	 0.0	 0.0	 1	 -30.0	 30.0;
	9	 10	 0.0	 0.11	 0.0	 267	 267	 267	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 12	 0.0	 0.256	 0.0	 115	 115	 115	 0.932	 0.0	 1	 -30.0	 30.0;
	12	 13	 0.0	 0.14	 0.0	 210	 210	 210	 0.0	 0.0	 1	 -30.0	 30.0;
	12	 14	 0.1231	 0.2559	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	12	 15	 0.0662	 0.1304	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	12	 16	 0.0945	 0.1987	 0.0	 30	 30	 30	 0.0	 0.0	 1	 -30.0	 30.0;
	14	 15	 0.221	 0.1997	 0.0	 20	 20	 20	 0.0	 0.0	 1	 -30.0	 30.0;
	16	 17	 0.0524	 0.1923	 0.0	 38	 38	 38	 0.0	 0.0	 1	 -30.0	 30.0;
	15	 18	 0.1073	 0.2185	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	18	 19	 0.0639	 0.1292	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	19	 20	 0.034	 0.068	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	10	 20	 0.0936	 0.209	 0.0	 30	 30	 30	 0.0	 0.0	 1	 -30.0	 30.0;
	10	 17	 0.0324	 0.0845	 0.0	 33	 33	 33	 0.0	 0.0	 1	 -30.0	 30.0;
	10	 21	 0.0348	 0.0749	 0.0	 30	 30	 30	 0.0	 0.0	 1	 -30.0	 30.0;
	10	 22	 0.0727	 0.1499	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	21	 22	 0.0116	 0.0236	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	15	 23	 0.1	 0.202	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	22	 24	 0.115	 0.179	 0.0	 26	 26	 26	 0.0	 0.0	 1	 -30.0	 30.0;
	23	 24	 0.132	 0.27	 0.0	 29	 29	 29	 0.0	 0.0	 1	 -30.0	 30.0;
	24	 25	 0.1885	 0.3292	 0.0	 27	 27	 27	 0.0	 0.0	 1	 -30.0	 30.0;
	25	 26	 0.2544	 0.38	 0.0	 25	 25	 25	 0.0	 0.0	 1	 -30.0	 30.0;
	25	 27	 0.1093	 0.2087	 0.0	 28	 28	 28	 0.0	 0.0	 1	 -30.0	 30.0;
	28	 27	 0.0	 0.396	 0.0	 75	 75	 75	 0.968	 0.0	 1	 -30.0	 30.0;
	27	 29	 0.2198	 0.4153	 0.0	 28	 28	 28	 0.0	 0.0	 1	 -30.0	 30.0;
	27	 30	 0.3202	 0.6027	 0.0	 28	 28	 28	 0.0	 0.0	 1	 -30.0	 30.0;
	29	 30	 0.2399	 0.4533	 0.0	 28	 28	 28	 0.0	 0.0	 1	 -30.0	 30.0;
	8	 28	 0.0636	 0.2	 0.0428	 140	 140	 140	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 28	 0.0169	 0.0599	 0.013	 149	 149	 149	 0.0	 0.0	 1	 -30.0	 30.0;
];


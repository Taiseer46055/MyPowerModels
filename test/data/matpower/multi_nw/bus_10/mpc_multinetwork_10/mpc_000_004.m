function mpc = mpc_000_004
%MPC_000_004	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 0 	Weight: 44
%	Time step: 4

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 44;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1695.6386492017637	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6586.39012375593	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5522.431315423298	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2618.601134544284	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	4523.335212213056	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	4509.081508269985	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3756.696946815666	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3926.734529003823	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	7007.577155349592	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4716.765915184095	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.46559093483552	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.506321816252624	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.13684185301526375	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.870032556344666	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0035058220491777045	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.572643439419882	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.04298775642064995	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.088620116406762	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.28888200994319	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.054204143464339535	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7288230828027005	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.042471023395659085	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.240399675428578	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.004648417090970509	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.109298196943538	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0010056709014922625	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.247561327758032	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.06119846280212126	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.797588944334441	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.10361090716794906	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.576120618375311	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.357153258075213	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.038028080736177215	0.0	0.0	221	221	2191	28.0	5;	%9 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1776.8	266.52	3.5	1	1	1	1776.8	7;	%DE1 0 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-266.52	-1776.8	3.5	1	1	1	-1776.8	8;	%DE1 0 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1294.0	194.1	3.5	1	1	1	1294.0	7;	%DE1 12 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-194.1	-1294.0	3.5	1	1	1	-1294.0	8;	%DE1 12 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.42	1.113	3.5	1	1	1	7.42	7;	%DE1 14 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-1.113	-7.42	3.5	1	1	1	-7.42	8;	%DE1 14 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.0	33.0	3.5	1	1	1	220.0	7;	%DE1 17 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-33.0	-220.0	3.5	1	1	1	-220.0	8;	%DE1 17 PHS (pump mode)
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	960.0	144.0	3.5	1	1	1	960.0	7;	%DE1 2 PHS
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	-144.0	-960.0	3.5	1	1	1	-960.0	8;	%DE1 2 PHS (pump mode)
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	39.8	5.97	3.5	1	1	1	39.8	7;	%DE1 20 PHS
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	-5.97	-39.8	3.5	1	1	1	-39.8	8;	%DE1 20 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	160.0	24.0	3.5	1	1	1	160.0	7;	%DE1 21 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.0	-160.0	3.5	1	1	1	-160.0	8;	%DE1 21 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	124.0	18.599999999999998	3.5	1	1	1	124.0	6;	%DE1 25 hydro
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	167.1	25.064999999999998	3.5	1	1	1	167.1	7;	%DE1 27 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-25.064999999999998	-167.1	3.5	1	1	1	-167.1	8;	%DE1 27 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	109.0	16.349999999999998	3.5	1	1	1	109.0	7;	%DE1 31 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-16.349999999999998	-109.0	3.5	1	1	1	-109.0	8;	%DE1 31 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	480.0	72.0	3.5	1	1	1	480.0	7;	%DE1 33 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-72.0	-480.0	3.5	1	1	1	-480.0	8;	%DE1 33 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	92.0	13.799999999999999	3.5	1	1	1	92.0	7;	%DE1 36 PHS
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.799999999999999	-92.0	3.5	1	1	1	-92.0	8;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	278.0	41.699999999999996	3.5	1	1	1	278.0	7;	%DE1 43 PHS
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	-41.699999999999996	-278.0	3.5	1	1	1	-278.0	8;	%DE1 43 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	164.0	24.599999999999998	3.5	1	1	1	164.0	7;	%DE1 49 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.599999999999998	-164.0	3.5	1	1	1	-164.0	8;	%DE1 49 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.7	11.955	3.5	1	1	1	79.7	7;	%DE1 50 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-11.955	-79.7	3.5	1	1	1	-79.7	8;	%DE1 50 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	399.8	59.97	3.5	1	1	1	399.8	7;	%DE1 52 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-59.97	-399.8	3.5	1	1	1	-399.8	8;	%DE1 52 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	162.0	24.3	3.5	1	1	1	162.0	7;	%DE1 54 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.3	-162.0	3.5	1	1	1	-162.0	8;	%DE1 54 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	90.0	13.5	3.5	1	1	1	90.0	7;	%DE1 61 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.5	-90.0	3.5	1	1	1	-90.0	8;	%DE1 61 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	1060.0	159.0	3.5	1	1	1	1060.0	7;	%DE1 7 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-159.0	-1060.0	3.5	1	1	1	-1060.0	8;	%DE1 7 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	45.5	6.825	3.5	1	1	1	45.5	6;	%DE1 71 hydro
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	18.0	3.5	1	1	1	120.0	7;	%DE1 76 PHS
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	-18.0	-120.0	3.5	1	1	1	-120.0	8;	%DE1 76 PHS (pump mode)
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0	investment
mpc.gencost = [
	2	366958.2302565932	146783.29210263726	3	0.0	1872.2358686560874	0.0	0.0	497.39955352724974;	%0 CCGT
	2	495032.9685009338	198013.18740037348	3	0.0	2525.678410719049	0.0	0.0	226.50678386591207;	%0 OCGT
	2	3350.332215788047	1340.132886315219	3	0.0	655.770643137218	0.0	0.0	1396.7642552130096;	%0 biomass
	2	0.0	0.0	3	0.0	1.087164073837639	0.0	0.0	966.4661572939276;	%0 offwind
	2	0.0	0.0	3	0.0	1.081398507661739	0.0	0.0	482.62318183544716;	%0 onwind
	2	0.0	0.0	3	0.0	0.865372720977405	0.0	0.0	178.82318889582143;	%0 solar
	2	368606.65230522206	147442.66092208883	3	0.0	1880.6461852307248	0.0	0.0	497.39955352724974;	%1 CCGT
	2	474102.2719742295	189640.9087896918	3	0.0	2418.8891427256603	0.0	0.0	226.50678386591207;	%1 OCGT
	2	3350.261664444093	1340.1046657776371	3	0.0	655.7568339095895	0.0	0.0	1396.7642552130096;	%1 biomass
	2	0.0	0.0	3	0.0	1.10064133779278	0.0	0.0	482.62318183544716;	%1 onwind
	2	0.0	0.0	3	0.0	0.885748615091563	0.0	0.0	178.82318889582143;	%1 solar
	2	351393.72623642907	140557.49049457163	3	0.0	1792.8251338593316	0.0	0.0	497.39955352724974;	%2 CCGT
	2	487134.5686357744	194853.8274543098	3	0.0	2485.3804522233386	0.0	0.0	226.50678386591207;	%2 OCGT
	2	3350.234038516853	1340.0936154067413	3	0.0	655.7514266034161	0.0	0.0	1396.7642552130092;	%2 biomass
	2	0.0	0.0	3	0.0	1.1049726749006477	0.0	0.0	482.62318183544716;	%2 onwind
	2	0.0	0.0	3	0.0	0.8782751721317004	0.0	0.0	178.82318889582143;	%2 solar
	2	325375.11795140884	130150.04718056355	3	0.0	1660.077132405147	0.0	0.0	497.39955352724974;	%3 CCGT
	2	467414.4398110065	186965.77592440258	3	0.0	2384.767550056155	0.0	0.0	226.50678386591207;	%3 OCGT
	2	3350.3294205900133	1340.1317682360052	3	0.0	655.7700960246649	0.0	0.0	1396.7642552130096;	%3 biomass
	2	0.0	0.0	3	0.0	1.1063742805806778	0.0	0.0	989.578376230859;	%3 offwind
	2	0.0	0.0	3	0.0	1.0865757159173526	0.0	0.0	482.62318183544727;	%3 onwind
	2	0.0	0.0	3	0.0	0.8734139436145897	0.0	0.0	178.82318889582143;	%3 solar
	2	348407.3533142363	139362.9413256945	3	0.0	1777.5885373175317	0.0	0.0	497.39955352724974;	%4 CCGT
	2	480855.6003100033	192342.2401240013	3	0.0	2453.344899540833	0.0	0.0	226.5067838659121;	%4 OCGT
	2	3350.333051537776	1340.1332206151103	3	0.0	655.7708067210365	0.0	0.0	1396.7642552130096;	%4 biomass
	2	0.0	0.0	3	0.0	1.0894746822841843	0.0	0.0	482.62318183544716;	%4 onwind
	2	0.0	0.0	3	0.0	0.8604376457232631	0.0	0.0	178.82318889582143;	%4 solar
	2	362754.98055548815	145101.99222219526	3	0.0	1850.7907171198374	0.0	0.0	497.39955352724974;	%5 CCGT
	2	483085.0589471374	193234.02357885498	3	0.0	2464.7196885058033	0.0	0.0	226.50678386591207;	%5 OCGT
	2	3350.3491813132932	1340.1396725253173	3	0.0	655.7739638507131	0.0	0.0	1396.7642552130096;	%5 biomass
	2	0.0	0.0	3	0.0	1.100775036031925	0.0	0.0	482.62318183544716;	%5 onwind
	2	0.0	0.0	3	0.0	0.8906291128484016	0.0	0.0	178.82318889582143;	%5 solar
	2	345029.37635347206	138011.75054138884	3	0.0	1760.3539609871023	0.0	0.0	497.39955352724974;	%6 CCGT
	2	487364.041120732	194945.6164482928	3	0.0	2486.551230207816	0.0	0.0	226.5067838659121;	%6 OCGT
	2	3350.3022630989867	1340.1209052395945	3	0.0	655.7647804069263	0.0	0.0	1396.7642552130096;	%6 biomass
	2	0.0	0.0	3	0.0	1.1162113250717252	0.0	0.0	482.62318183544716;	%6 onwind
	2	0.0	0.0	3	0.0	0.8888750302277983	0.0	0.0	178.82318889582137;	%6 solar
	2	350644.4614070975	140257.78456283902	3	0.0	1789.0023541178443	0.0	0.0	497.39955352724974;	%7 CCGT
	2	490615.8735558812	196246.34942235248	3	0.0	2503.1422120198017	0.0	0.0	226.50678386591207;	%7 OCGT
	2	3350.346174258694	1340.1384697034775	3	0.0	655.7733752708345	0.0	0.0	1396.7642552130096;	%7 biomass
	2	0.0	0.0	3	0.0	1.107576004574105	0.0	0.0	482.62318183544716;	%7 onwind
	2	0.0	0.0	3	0.0	0.8732419688836109	0.0	0.0	178.82318889582143;	%7 solar
	2	377113.85587496654	150845.54234998664	3	0.0	1924.05028507636	0.0	0.0	497.39955352724974;	%8 CCGT
	2	493314.3589007459	197325.74356029835	3	0.0	2516.9099943915603	0.0	0.0	226.50678386591207;	%8 OCGT
	2	3350.372288417991	1340.1489153671964	3	0.0	655.7784866741027	0.0	0.0	1396.7642552130096;	%8 biomass
	2	0.0	0.0	3	0.0	1.115032523322936	0.0	0.0	482.62318183544716;	%8 onwind
	2	0.0	0.0	3	0.0	0.8777811093152537	0.0	0.0	178.82318889582143;	%8 solar
	2	379780.5163080527	151912.20652322107	3	0.0	1937.6556954492482	0.0	0.0	497.39955352724974;	%9 CCGT
	2	495171.00268819893	198068.4010752796	3	0.0	2526.382666776525	0.0	0.0	226.5067838659121;	%9 OCGT
	2	3350.240433693778	1340.096173477511	3	0.0	655.7526783507101	0.0	0.0	1396.7642552130096;	%9 biomass
	2	0.0	0.0	3	0.0	1.0900695123226547	0.0	0.0	1022.0314052199734;	%9 offwind
	2	0.0	0.0	3	0.0	1.0990467253515237	0.0	0.0	482.62318183544716;	%9 onwind
	2	0.0	0.0	3	0.0	0.879621122472658	0.0	0.0	178.82318889582143;	%9 solar
	2	0.0	0.0	3	0.0	0.42539582476808047	0.0	0.0	7803189.531257401;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.31904686857606035	0.0	0.0	7803189.531257401;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4634606411349178	0.0	0.0	7803189.531257401;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.34759548085118835	0.0	0.0	7803189.531257401;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.48359632668291186	0.0	0.0	7803189.531257401;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.3626972450121839	0.0	0.0	7803189.531257401;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4038515402033509	0.0	0.0	7803189.531257401;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.3028886551525132	0.0	0.0	7803189.531257401;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4804662147609432	0.0	0.0	7803189.531257401;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.3603496610707074	0.0	0.0	7803189.531257401;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4622088238105458	0.0	0.0	7803189.531257401;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.34665661785790935	0.0	0.0	7803189.531257401;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4399614464380775	0.0	0.0	7803189.531257401;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.32997108482855814	0.0	0.0	7803189.531257401;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.44659383921462215	0.0	0.0	7803189.531257401;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.44156727639977206	0.0	0.0	7803189.531257401;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.33117545729982906	0.0	0.0	7803189.531257401;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4299871542278731	0.0	0.0	7803189.531257401;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.32249036567090483	0.0	0.0	7803189.531257401;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.422083852074183	0.0	0.0	7803189.531257401;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.31656288905563723	0.0	0.0	7803189.531257401;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40420711277944493	0.0	0.0	7803189.531257401;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.3031553345845837	0.0	0.0	7803189.531257401;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.44582789170418347	0.0	0.0	7803189.531257401;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.3343709187781376	0.0	0.0	7803189.531257401;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40924822064297683	0.0	0.0	7803189.531257401;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.3069361654822326	0.0	0.0	7803189.531257401;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4117103418935138	0.0	0.0	7803189.531257401;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.30878275642013536	0.0	0.0	7803189.531257401;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.47711994389072127	0.0	0.0	7803189.531257401;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.35783995791804096	0.0	0.0	7803189.531257401;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43676627730082007	0.0	0.0	7803189.531257401;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.32757470797561505	0.0	0.0	7803189.531257401;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4539344736067687	0.0	0.0	7803189.531257401;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.3404508552050765	0.0	0.0	7803189.531257401;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.460197860095951	0.0	0.0	7803189.531257401;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.34514839507196327	0.0	0.0	7803189.531257401;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40860130559967683	0.0	0.0	7803189.531257401;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.475607582503133	0.0	0.0	7803189.531257401;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.35670568687734977	0.0	0.0	7803189.531257401;	%DE1 76 PHS (pump mode)
];
%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0.002346193481236678	0.019238786546140762	1.958387880013024e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	1	8	0.0023321905176023603	0.019123962244339356	3.236859740616215e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	1	9	0.0019445334507215365	0.0159451742959166	2.6988284162105606e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	2	3	0.002476449448812836	0.020306885480265262	2.067113656570062e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	2	6	0.00032235809806776487	0.002643336404155672	9.943293073955649e-05	20645.352805897994	20645.352805897994	20645.352805897994	0	0	1	-60	60;
	2	7	0.00018528201517929854	0.001519312524470248	8.729978747518073e-05	25516.226086943192	25516.226086943192	25516.226086943192	0	0	1	-60	60;
	2	8	0.0027328436533095644	0.02240931795713843	2.2811281044864822e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	2	10	0.000789978277571624	0.006477821876087317	8.077670377089266e-05	11886.718282183694	11886.718282183694	11886.718282183694	0	0	1	-60	60;
	3	5	0.0007685374454869242	0.006302007052992778	4.9916538212103785e-05	9473.625097078733	9473.625097078733	9473.625097078733	0	0	1	-60	60;
	3	6	0.0011358075415830791	0.00931362184098125	4.969480021174222e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	3	8	0.0007808106224918078	0.006402647104432824	3.4162678506126354e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	3	9	0.0009715277760101635	0.007966527763283342	4.250709469818385e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	4	8	0.002553846799021281	0.0209415437519745	2.1317178905369238e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	4	10	0.0007315818716301728	0.0059989713473674175	4.061465629624203e-05	8758.634523714301	8758.634523714301	8758.634523714301	0	0	1	-60	60;
	5	6	0.0007449955480032163	0.0061089634936263735	4.838749101703001e-05	9473.625097078733	9473.625097078733	9473.625097078733	0	0	1	-60	60;
	7	10	0.00047207592830423385	0.0038710226120947182	6.304733408087732e-05	13584.820893924221	13584.820893924221	13584.820893924221	0	0	1	-60	60;
	8	9	0.0005598156111912242	0.00459048801176804	5.056298758014288e-05	11171.727708819262	11171.727708819262	11171.727708819262	0	0	1	-60	60;
	8	10	0.0008717140618085548	0.0071480553068301485	3.813996170774132e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
];

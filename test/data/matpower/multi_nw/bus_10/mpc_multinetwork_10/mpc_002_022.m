function mpc = mpc_002_022
%MPC_002_022	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 2 	Weight: 96
%	Time step: 22

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 96;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	339.4761059766178	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5114.511425221561	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4187.307198955227	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	1651.0913559173187	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2187.658121178794	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-133.94135496903968	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8298.136969445508	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3805.1001096161053	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4266.916120570282	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	6496.317022078017	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511627	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.6919883178594395	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.08273272111183	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7786889635120304	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.3666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.346391820006516	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.176840437460433	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.013686465360125	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.55499999999998	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.011805068448621	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.80000000000004	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.254524630441825	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49474999999998	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.229187992126858	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.65999999999997	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781609	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.556737346345655	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.350278072309181	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.0553846153846	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5774802428990324	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.073443918492178	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	800636.1387416578	320254.4554966631	3	0.0	4084.878258886009	0.0	0.0	1085.2353895139995;	%0 CCGT
	2	1080071.9312747647	432028.7725099058	3	0.0	5510.571077932471	0.0	0.0	494.19661934380815;	%0 OCGT
	2	7309.815743537557	2923.926297415023	3	0.0	1430.7723122993848	0.0	0.0	3047.485647737475;	%0 biomass
	2	0.0	0.0	3	0.0	2.371994342918485	0.0	0.0	2108.653434095842;	%0 offwind
	2	0.0	0.0	3	0.0	2.359414925807431	0.0	0.0	1052.9960330955212;	%0 onwind
	2	0.0	0.0	3	0.0	1.8880859366779745	0.0	0.0	390.1596848636104;	%0 solar
	2	804232.6959386664	321693.0783754665	3	0.0	4103.2280405034	0.0	0.0	1085.2353895139995;	%1 CCGT
	2	1034404.9570346826	413761.982813873	3	0.0	5277.576311401441	0.0	0.0	494.19661934380815;	%1 OCGT
	2	7309.661813332566	2923.8647253330264	3	0.0	1430.742183075468	0.0	0.0	3047.485647737475;	%1 biomass
	2	0.0	0.0	3	0.0	2.4013992824569743	0.0	0.0	1052.9960330955212;	%1 onwind
	2	0.0	0.0	3	0.0	1.9325424329270464	0.0	0.0	390.1596848636104;	%1 solar
	2	766677.2208794816	306670.8883517926	3	0.0	3911.6184738749052	0.0	0.0	1085.2353895139995;	%2 CCGT
	2	1062839.0588416897	425135.62353667594	3	0.0	5422.648259396376	0.0	0.0	494.19661934380815;	%2 OCGT
	2	7309.601538582225	2923.8406154328904	3	0.0	1430.7303853165445	0.0	0.0	3047.4856477374747;	%2 biomass
	2	0.0	0.0	3	0.0	2.410849472510504	0.0	0.0	1052.9960330955212;	%2 onwind
	2	0.0	0.0	3	0.0	1.9162367391964372	0.0	0.0	390.1596848636104;	%2 solar
	2	709909.3482576193	283963.7393030477	3	0.0	3621.986470702139	0.0	0.0	1085.2353895139995;	%3 CCGT
	2	1019813.3232240141	407925.3292896056	3	0.0	5203.12920012252	0.0	0.0	494.19661934380815;	%3 OCGT
	2	7309.8096449236655	2923.923857969466	3	0.0	1430.771118599269	0.0	0.0	3047.485647737475;	%3 biomass
	2	0.0	0.0	3	0.0	2.4139075212669336	0.0	0.0	2159.0800935946013;	%3 offwind
	2	0.0	0.0	3	0.0	2.3707106529105872	0.0	0.0	1052.9960330955214;	%3 onwind
	2	0.0	0.0	3	0.0	1.9056304224318321	0.0	0.0	390.1596848636104;	%3 solar
	2	760161.4981401518	304064.5992560607	3	0.0	3878.3749905109785	0.0	0.0	1085.2353895139995;	%4 CCGT
	2	1049139.491585462	419655.79663418466	3	0.0	5352.75250808909	0.0	0.0	494.19661934380827;	%4 OCGT
	2	7309.817566991512	2923.9270267966044	3	0.0	1430.7726692095343	0.0	0.0	3047.485647737475;	%4 biomass
	2	0.0	0.0	3	0.0	2.3770356704382203	0.0	0.0	1052.9960330955212;	%4 onwind
	2	0.0	0.0	3	0.0	1.8773184997598467	0.0	0.0	390.1596848636104;	%4 solar
	2	791465.4121210651	316586.164848426	3	0.0	4038.0888373523726	0.0	0.0	1085.2353895139995;	%5 CCGT
	2	1054003.7649755725	421601.505990229	3	0.0	5377.570229467206	0.0	0.0	494.19661934380815;	%5 OCGT
	2	7309.852759229003	2923.9411036916013	3	0.0	1430.7795574924648	0.0	0.0	3047.485647737475;	%5 biomass
	2	0.0	0.0	3	0.0	2.4016909877060186	0.0	0.0	1052.9960330955212;	%5 onwind
	2	0.0	0.0	3	0.0	1.94319079166924	0.0	0.0	390.1596848636104;	%5 solar
	2	752791.3665893937	301116.54663575743	3	0.0	3840.772278517314	0.0	0.0	1085.2353895139995;	%6 CCGT
	2	1063339.7260815972	425335.8904326388	3	0.0	5425.202684089781	0.0	0.0	494.19661934380827;	%6 OCGT
	2	7309.750392215971	2923.9001568863882	3	0.0	1430.7595208878392	0.0	0.0	3047.485647737475;	%6 biomass
	2	0.0	0.0	3	0.0	2.435370163792855	0.0	0.0	1052.9960330955212;	%6 onwind
	2	0.0	0.0	3	0.0	1.9393637023151964	0.0	0.0	390.1596848636103;	%6 solar
	2	765042.4612518492	306016.98450073967	3	0.0	3903.2778635298423	0.0	0.0	1085.2353895139995;	%7 CCGT
	2	1070434.6332128318	428173.8532851327	3	0.0	5461.401189861386	0.0	0.0	494.19661934380815;	%7 OCGT
	2	7309.8461983826055	2923.938479353042	3	0.0	1430.7782733181846	0.0	0.0	3047.485647737475;	%7 biomass
	2	0.0	0.0	3	0.0	2.41652946452532	0.0	0.0	1052.9960330955212;	%7 onwind
	2	0.0	0.0	3	0.0	1.9052552048369695	0.0	0.0	390.1596848636104;	%7 solar
	2	822793.8673635634	329117.5469454254	3	0.0	4197.927894712058	0.0	0.0	1085.2353895139995;	%8 CCGT
	2	1076322.2376016274	430528.89504065097	3	0.0	5491.439987763404	0.0	0.0	494.19661934380815;	%8 OCGT
	2	7309.9031747301615	2923.9612698920646	3	0.0	1430.7894254707694	0.0	0.0	3047.485647737475;	%8 biomass
	2	0.0	0.0	3	0.0	2.4327982327045876	0.0	0.0	1052.9960330955212;	%8 onwind
	2	0.0	0.0	3	0.0	1.9151587839605535	0.0	0.0	390.1596848636104;	%8 solar
	2	828612.0355812059	331444.81423248234	3	0.0	4227.612426434724	0.0	0.0	1085.2353895139995;	%9 CCGT
	2	1080373.0967742521	432149.23870970093	3	0.0	5512.107636603328	0.0	0.0	494.19661934380827;	%9 OCGT
	2	7309.615491695515	2923.846196678206	3	0.0	1430.733116401549	0.0	0.0	3047.485647737475;	%9 biomass
	2	0.0	0.0	3	0.0	2.3783334814312465	0.0	0.0	2229.886702298124;	%9 offwind
	2	0.0	0.0	3	0.0	2.3979201280396882	0.0	0.0	1052.9960330955212;	%9 onwind
	2	0.0	0.0	3	0.0	1.9191733581221628	0.0	0.0	390.1596848636104;	%9 solar
	2	0.0	0.0	3	0.0	0.9281363449485391	0.0	0.0	17025140.795470692;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.6961022587114043	0.0	0.0	17025140.795470692;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.0111868533852753	0.0	0.0	17025140.795470692;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.7583901400389564	0.0	0.0	17025140.795470692;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.0551192582172622	0.0	0.0	17025140.795470692;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.7913394436629466	0.0	0.0	17025140.795470692;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.8811306331709474	0.0	0.0	17025140.795470692;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.6608479748782106	0.0	0.0	17025140.795470692;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.048289923114785	0.0	0.0	17025140.795470692;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.7862174423360888	0.0	0.0	17025140.795470692;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.0084556155866453	0.0	0.0	17025140.795470692;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.756341711689984	0.0	0.0	17025140.795470692;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9599158831376237	0.0	0.0	17025140.795470692;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.7199369123532178	0.0	0.0	17025140.795470692;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9743865582864484	0.0	0.0	17025140.795470692;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.9634195121449572	0.0	0.0	17025140.795470692;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.722564634108718	0.0	0.0	17025140.795470692;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9381537910426323	0.0	0.0	17025140.795470692;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.7036153432819743	0.0	0.0	17025140.795470692;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9209102227073083	0.0	0.0	17025140.795470692;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.6906826670304812	0.0	0.0	17025140.795470692;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.8819064278824253	0.0	0.0	17025140.795470692;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.661429820911819	0.0	0.0	17025140.795470692;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9727154000818549	0.0	0.0	17025140.795470692;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.7295365500613912	0.0	0.0	17025140.795470692;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.8929052086755858	0.0	0.0	17025140.795470692;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.6696789065066893	0.0	0.0	17025140.795470692;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.8982771095858484	0.0	0.0	17025140.795470692;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.6737078321893863	0.0	0.0	17025140.795470692;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.0409889684888465	0.0	0.0	17025140.795470692;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.7807417263666349	0.0	0.0	17025140.795470692;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.952944605019971	0.0	0.0	17025140.795470692;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.7147084537649783	0.0	0.0	17025140.795470692;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.9904024878693136	0.0	0.0	17025140.795470692;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.7428018659019853	0.0	0.0	17025140.795470692;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	1.0040680583911659	0.0	0.0	17025140.795470692;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.7530510437933744	0.0	0.0	17025140.795470692;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.8914937576720222	0.0	0.0	17025140.795470692;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	1.0376892709159267	0.0	0.0	17025140.795470692;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.778266953186945	0.0	0.0	17025140.795470692;	%DE1 76 PHS (pump mode)
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

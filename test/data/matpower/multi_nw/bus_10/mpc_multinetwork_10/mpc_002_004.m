function mpc = mpc_002_004
%MPC_002_004	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 2 	Weight: 21
%	Time step: 4

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 21;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1172.342901409961	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6413.247933382012	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5135.165348366134	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2400.809305601276	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	4259.924201325188	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	2560.5585377133007	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3930.1620326930583	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3705.6110925860303	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	7166.841177861013	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4246.486795714735	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.471748459219693	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.295354471629018	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.286735529474786	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.638107371319781	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.032363249371317	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.067428625924245	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8058257162101226	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6601823267027442	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.023848882573482	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.45206478584889	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.52524375395798	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.71070025639359	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.932314309341715	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	175139.15534973764	70055.66213989505	3	0.0	893.5671191313145	0.0	0.0	53176.53408618597;	%0 CCGT
	2	236265.73496635476	94506.29398654189	3	0.0	1205.437423297728	0.0	0.0	24215.634347846604;	%0 OCGT
	2	1599.0221938988407	639.6088775595363	3	0.0	312.9814433154904	0.0	0.0	6666.374854425727;	%0 biomass
	2	0.0	0.0	3	0.0	0.5188737625134185	0.0	0.0	7841.554958043913;	%0 offwind
	2	0.0	0.0	3	0.0	0.5161220150203755	0.0	0.0	8753.029525106518;	%0 onwind
	2	0.0	0.0	3	0.0	0.41301879864830693	0.0	0.0	2389.7280697896135;	%0 solar
	2	175925.90223658326	70370.3608946333	3	0.0	897.5811338601186	0.0	0.0	53176.53408618597;	%1 CCGT
	2	226276.08435133682	90510.43374053472	3	0.0	1154.4698181190652	0.0	0.0	24215.634347846604;	%1 OCGT
	2	1598.9885216664989	639.5954086665995	3	0.0	312.97485254775864	0.0	0.0	6666.374854425727;	%1 biomass
	2	0.0	0.0	3	0.0	0.5253060930374631	0.0	0.0	8753.029525106518;	%1 onwind
	2	0.0	0.0	3	0.0	0.4227436572027914	0.0	0.0	2389.7280697896135;	%1 solar
	2	167710.6420673866	67084.25682695463	3	0.0	855.6665411601355	0.0	0.0	53176.53408618597;	%2 CCGT
	2	232496.04412161963	92998.41764864785	3	0.0	1186.2043067429572	0.0	0.0	24215.634347846604;	%2 OCGT
	2	1598.9753365648617	639.5901346259448	3	0.0	312.9722717879941	0.0	0.0	6666.374854425727;	%2 biomass
	2	0.0	0.0	3	0.0	0.5273733221116728	0.0	0.0	8753.029525106518;	%2 onwind
	2	0.0	0.0	3	0.0	0.41917678669922065	0.0	0.0	2389.7280697896135;	%2 solar
	2	155292.66993135423	62117.06797254169	3	0.0	792.3095404660929	0.0	0.0	53176.53408618597;	%3 CCGT
	2	223084.1644552531	89233.66578210123	3	0.0	1138.1845125268012	0.0	0.0	24215.634347846604;	%3 OCGT
	2	1599.020859827052	639.6083439308206	3	0.0	312.9811821935901	0.0	0.0	6666.374854425727;	%3 biomass
	2	0.0	0.0	3	0.0	0.5280422702771417	0.0	0.0	8029.079098054924;	%3 offwind
	2	0.0	0.0	3	0.0	0.518592955324191	0.0	0.0	8753.029525106518;	%3 onwind
	2	0.0	0.0	3	0.0	0.4168566549069633	0.0	0.0	2389.7280697896135;	%3 solar
	2	166285.32771815822	66514.13108726329	3	0.0	848.3945291742765	0.0	0.0	53176.53408618597;	%4 CCGT
	2	229499.26378431977	91799.7055137279	3	0.0	1170.9146111444884	0.0	0.0	24215.634347846604;	%4 OCGT
	2	1599.022592779393	639.6090371117572	3	0.0	312.9815213895856	0.0	0.0	6666.374854425727;	%4 biomass
	2	0.0	0.0	3	0.0	0.5199765529083606	0.0	0.0	8753.029525106518;	%4 onwind
	2	0.0	0.0	3	0.0	0.41066342182246646	0.0	0.0	2389.7280697896135;	%4 solar
	2	173133.058901483	69253.2235605932	3	0.0	883.3319331708315	0.0	0.0	53176.53408618597;	%5 CCGT
	2	230563.3235884065	92225.3294353626	3	0.0	1176.3434876959516	0.0	0.0	24215.634347846604;	%5 OCGT
	2	1599.0302910813443	639.6121164325378	3	0.0	312.9830282014767	0.0	0.0	6666.374854425727;	%5 biomass
	2	0.0	0.0	3	0.0	0.5253699035606916	0.0	0.0	8753.029525106518;	%5 onwind
	2	0.0	0.0	3	0.0	0.42507298567764623	0.0	0.0	2389.7280697896135;	%5 solar
	2	164673.11144142985	65869.24457657195	3	0.0	840.1689359256625	0.0	0.0	53176.53408618597;	%6 CCGT
	2	232605.5650803494	93042.22603213975	3	0.0	1186.7630871446395	0.0	0.0	24215.634347846604;	%6 OCGT
	2	1599.0078982972436	639.6031593188974	3	0.0	312.97864519421483	0.0	0.0	6666.374854425727;	%6 biomass
	2	0.0	0.0	3	0.0	0.532737223329687	0.0	0.0	8753.029525106518;	%6 onwind
	2	0.0	0.0	3	0.0	0.4242358098814492	0.0	0.0	2389.7280697896135;	%6 solar
	2	167353.038398842	66941.2153595368	3	0.0	853.842032647153	0.0	0.0	53176.53408618597;	%7 CCGT
	2	234157.57601530693	93663.03040612277	3	0.0	1194.6815102821781	0.0	0.0	24215.634347846604;	%7 OCGT
	2	1599.028855896195	639.6115423584779	3	0.0	312.98274728835287	0.0	0.0	6666.374854425727;	%7 biomass
	2	0.0	0.0	3	0.0	0.5286158203649137	0.0	0.0	8753.029525106518;	%7 onwind
	2	0.0	0.0	3	0.0	0.41677457605808704	0.0	0.0	2389.7280697896135;	%7 solar
	2	179986.15848577948	71994.46339431181	3	0.0	918.2967269682628	0.0	0.0	53176.53408618597;	%8 CCGT
	2	235445.489475356	94178.1957901424	3	0.0	1201.2524973232446	0.0	0.0	24215.634347846604;	%8 OCGT
	2	1599.041319472223	639.6165277888891	3	0.0	312.98518682173085	0.0	0.0	6666.374854425727;	%8 biomass
	2	0.0	0.0	3	0.0	0.5321746134041285	0.0	0.0	8753.029525106518;	%8 onwind
	2	0.0	0.0	3	0.0	0.4189409839913711	0.0	0.0	2389.7280697896135;	%8 solar
	2	181258.88278338878	72503.5531133555	3	0.0	924.7902182825958	0.0	0.0	53176.53408618597;	%9 CCGT
	2	236331.61491936768	94532.64596774708	3	0.0	1205.7735455069778	0.0	0.0	24215.634347846604;	%9 OCGT
	2	1598.978388808394	639.5913555233576	3	0.0	312.9728692128389	0.0	0.0	6666.374854425727;	%9 biomass
	2	0.0	0.0	3	0.0	0.5202604490630852	0.0	0.0	8292.391174171147;	%9 offwind
	2	0.0	0.0	3	0.0	0.5245450280086817	0.0	0.0	8753.029525106518;	%9 onwind
	2	0.0	0.0	3	0.0	0.41981917208922315	0.0	0.0	2389.7280697896135;	%9 solar
	2	0.0	0.0	3	0.0	0.20302982545749293	0.0	0.0	425.14264258096057;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.1522723690931197	0.0	0.0	425.14264258096057;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22119712417802895	0.0	0.0	425.14264258096057;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.1658978431335217	0.0	0.0	425.14264258096057;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23080733773502612	0.0	0.0	425.14264258096057;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.1731055033012696	0.0	0.0	425.14264258096057;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19274732600614475	0.0	0.0	425.14264258096057;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.14456049450460856	0.0	0.0	425.14264258096057;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22931342068135926	0.0	0.0	425.14264258096057;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.17198506551101944	0.0	0.0	425.14264258096057;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22059966590957866	0.0	0.0	425.14264258096057;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.165449749432184	0.0	0.0	425.14264258096057;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20998159943635517	0.0	0.0	425.14264258096057;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.15748619957726637	0.0	0.0	425.14264258096057;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2131470596251606	0.0	0.0	425.14264258096057;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.2107480182817094	0.0	0.0	425.14264258096057;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.15806101371128203	0.0	0.0	425.14264258096057;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2052211417905758	0.0	0.0	425.14264258096057;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.15391585634293184	0.0	0.0	425.14264258096057;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2014491112172237	0.0	0.0	425.14264258096057;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.1510868334129178	0.0	0.0	425.14264258096057;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19291703109928054	0.0	0.0	425.14264258096057;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.14468777332446042	0.0	0.0	425.14264258096057;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21278149376790575	0.0	0.0	425.14264258096057;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.1595861203259293	0.0	0.0	425.14264258096057;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19532301439778438	0.0	0.0	425.14264258096057;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.14649226079833827	0.0	0.0	425.14264258096057;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19649811772190431	0.0	0.0	425.14264258096057;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.14737358829142824	0.0	0.0	425.14264258096057;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22771633685693515	0.0	0.0	425.14264258096057;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.17078725264270136	0.0	0.0	425.14264258096057;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20845663234811868	0.0	0.0	425.14264258096057;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.156342474261089	0.0	0.0	425.14264258096057;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21665054422141233	0.0	0.0	425.14264258096057;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.16248790816605924	0.0	0.0	425.14264258096057;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2196398877730675	0.0	0.0	425.14264258096057;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.16472991582980062	0.0	0.0	425.14264258096057;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19501425949075485	0.0	0.0	425.14264258096057;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.22699452801285894	0.0	0.0	425.14264258096057;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.1702458960096442	0.0	0.0	425.14264258096057;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_004_021
%MPC_004_021	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 4 	Weight: 27
%	Time step: 21

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 27;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1196.6084131632279	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6585.306545877138	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5217.428095021869	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2425.245256528967	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3499.0535981366916	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	3784.429558162405	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	1939.3745342078369	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3863.970662496386	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	7753.486234648054	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	2755.5580204785165	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.34897503406995	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.489232624716651	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.38385303527386	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.435324541963236	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.110228785469834	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.290481794511674	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.193625306575697	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.469167843492919	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.884728752658734	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.468180157718979	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.011390298645581	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.860657259105858	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.200407517276096	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	225178.91402109127	90071.5656084365	3	0.0	1148.87201031169	0.0	0.0	68369.82953938197;	%0 CCGT
	2	303770.23067102756	121508.09226841101	3	0.0	1549.8481156685075	0.0	0.0	31134.38701865992;	%0 OCGT
	2	2055.885677869938	822.3542711479752	3	0.0	402.404712834202	0.0	0.0	8571.053384261648;	%0 biomass
	2	0.0	0.0	3	0.0	0.6671234089458239	0.0	0.0	10081.999231770746;	%0 offwind
	2	0.0	0.0	3	0.0	0.6635854478833398	0.0	0.0	11253.89510370838;	%0 onwind
	2	0.0	0.0	3	0.0	0.5310241696906803	0.0	0.0	3072.507518300932;	%0 solar
	2	226190.44573274988	90476.17829309996	3	0.0	1154.032886391581	0.0	0.0	68369.82953938197;	%1 CCGT
	2	290926.3941660045	116370.55766640179	3	0.0	1484.3183375816552	0.0	0.0	31134.38701865992;	%1 OCGT
	2	2055.8423849997844	822.3369539999137	3	0.0	402.39623898997536	0.0	0.0	8571.053384261648;	%1 biomass
	2	0.0	0.0	3	0.0	0.675393548191024	0.0	0.0	11253.89510370838;	%1 onwind
	2	0.0	0.0	3	0.0	0.5435275592607318	0.0	0.0	3072.507518300932;	%1 solar
	2	215627.9683723542	86251.18734894168	3	0.0	1100.1426957773172	0.0	0.0	68369.82953938197;	%2 CCGT
	2	298923.48529922526	119569.3941196901	3	0.0	1525.1198229552306	0.0	0.0	31134.38701865992;	%2 OCGT
	2	2055.8254327262507	822.3301730905004	3	0.0	402.3929208702781	0.0	0.0	8571.053384261648;	%2 biomass
	2	0.0	0.0	3	0.0	0.6780514141435792	0.0	0.0	11253.89510370838;	%2 onwind
	2	0.0	0.0	3	0.0	0.538941582898998	0.0	0.0	3072.507518300932;	%2 solar
	2	199662.00419745545	79864.80167898217	3	0.0	1018.6836948849766	0.0	0.0	68369.82953938197;	%3 CCGT
	2	286822.49715675396	114728.99886270158	3	0.0	1463.3800875344589	0.0	0.0	31134.38701865992;	%3 OCGT
	2	2055.883962634781	822.3535850539123	3	0.0	402.40437710604436	0.0	0.0	8571.053384261648;	%3 biomass
	2	0.0	0.0	3	0.0	0.678911490356325	0.0	0.0	10323.101697499187;	%3 offwind
	2	0.0	0.0	3	0.0	0.6667623711311027	0.0	0.0	11253.895103708383;	%3 onwind
	2	0.0	0.0	3	0.0	0.5359585563089527	0.0	0.0	3072.507518300932;	%3 solar
	2	213795.42135191773	85518.16854076709	3	0.0	1090.7929660812126	0.0	0.0	68369.82953938197;	%4 CCGT
	2	295070.48200841114	118028.19280336445	3	0.0	1505.4616429000566	0.0	0.0	31134.38701865992;	%4 OCGT
	2	2055.8861907163628	822.354476286545	3	0.0	402.4048132151815	0.0	0.0	8571.053384261648;	%4 biomass
	2	0.0	0.0	3	0.0	0.6685412823107494	0.0	0.0	11253.89510370838;	%4 onwind
	2	0.0	0.0	3	0.0	0.5279958280574569	0.0	0.0	3072.507518300932;	%4 solar
	2	222599.64715904955	89039.85886361983	3	0.0	1135.7124855053548	0.0	0.0	68369.82953938197;	%5 CCGT
	2	296438.55889937974	118575.42355975191	3	0.0	1512.4416270376519	0.0	0.0	31134.38701865992;	%5 OCGT
	2	2055.896088533157	822.3584354132629	3	0.0	402.40675054475577	0.0	0.0	8571.053384261648;	%5 biomass
	2	0.0	0.0	3	0.0	0.6754755902923177	0.0	0.0	11253.89510370838;	%5 onwind
	2	0.0	0.0	3	0.0	0.5465224101569737	0.0	0.0	3072.507518300932;	%5 solar
	2	211722.57185326697	84689.02874130679	3	0.0	1080.2172033329946	0.0	0.0	68369.82953938197;	%6 CCGT
	2	299064.2979604492	119625.71918417967	3	0.0	1525.8382549002508	0.0	0.0	31134.38701865992;	%6 OCGT
	2	2055.867297810742	822.3469191242967	3	0.0	402.4011152497048	0.0	0.0	8571.053384261648;	%6 biomass
	2	0.0	0.0	3	0.0	0.6849478585667405	0.0	0.0	11253.89510370838;	%6 onwind
	2	0.0	0.0	3	0.0	0.545446041276149	0.0	0.0	3072.5075183009317;	%6 solar
	2	215168.1922270826	86067.27689083303	3	0.0	1097.7968991177681	0.0	0.0	68369.82953938197;	%7 CCGT
	2	301059.7405911089	120423.89623644357	3	0.0	1536.0190846485148	0.0	0.0	31134.38701865992;	%7 OCGT
	2	2055.894243295108	822.357697318043	3	0.0	402.4063893707394	0.0	0.0	8571.053384261648;	%7 biomass
	2	0.0	0.0	3	0.0	0.6796489118977462	0.0	0.0	11253.89510370838;	%7 onwind
	2	0.0	0.0	3	0.0	0.5358530263603977	0.0	0.0	3072.507518300932;	%7 solar
	2	231410.7751960022	92564.31007840089	3	0.0	1180.6672203877665	0.0	0.0	68369.82953938197;	%8 CCGT
	2	302715.62932545773	121086.25173018308	3	0.0	1544.4674965584575	0.0	0.0	31134.38701865992;	%8 OCGT
	2	2055.910267892858	822.3641071571432	3	0.0	402.4095259136539	0.0	0.0	8571.053384261648;	%8 biomass
	2	0.0	0.0	3	0.0	0.6842245029481653	0.0	0.0	11253.89510370838;	%8 onwind
	2	0.0	0.0	3	0.0	0.5386384079889057	0.0	0.0	3072.507518300932;	%8 solar
	2	233047.13500721415	93218.85400288567	3	0.0	1189.015994934766	0.0	0.0	68369.82953938197;	%9 CCGT
	2	303854.93346775847	121541.97338710338	3	0.0	1550.2802727946857	0.0	0.0	31134.38701865992;	%9 OCGT
	2	2055.8293570393635	822.3317428157454	3	0.0	402.3936889879357	0.0	0.0	8571.053384261648;	%9 biomass
	2	0.0	0.0	3	0.0	0.6689062916525381	0.0	0.0	10661.645795362903;	%9 offwind
	2	0.0	0.0	3	0.0	0.6744150360111623	0.0	0.0	11253.89510370838;	%9 onwind
	2	0.0	0.0	3	0.0	0.5397675069718584	0.0	0.0	3072.507518300932;	%9 solar
	2	0.0	0.0	3	0.0	0.26103834701677664	0.0	0.0	546.6119690326636;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.19577876026258248	0.0	0.0	546.6119690326636;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2843963025146086	0.0	0.0	546.6119690326636;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.21329722688595648	0.0	0.0	546.6119690326636;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29675229137360504	0.0	0.0	546.6119690326636;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.22256421853020378	0.0	0.0	546.6119690326636;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24781799057932896	0.0	0.0	546.6119690326636;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.18586349293449672	0.0	0.0	546.6119690326636;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29483154087603336	0.0	0.0	546.6119690326636;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.22112365565702502	0.0	0.0	546.6119690326636;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.283628141883744	0.0	0.0	546.6119690326636;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.212721106412808	0.0	0.0	546.6119690326636;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26997634213245664	0.0	0.0	546.6119690326636;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.2024822565993425	0.0	0.0	546.6119690326636;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2740462195180636	0.0	0.0	546.6119690326636;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.2709617377907692	0.0	0.0	546.6119690326636;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.20322130334307692	0.0	0.0	546.6119690326636;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26385575373074033	0.0	0.0	546.6119690326636;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.19789181529805525	0.0	0.0	546.6119690326636;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2590060001364305	0.0	0.0	546.6119690326636;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.19425450010232287	0.0	0.0	546.6119690326636;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2480361828419321	0.0	0.0	546.6119690326636;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.18602713713144908	0.0	0.0	546.6119690326636;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2735762062730217	0.0	0.0	546.6119690326636;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.20518215470476628	0.0	0.0	546.6119690326636;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2511295899400085	0.0	0.0	546.6119690326636;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.18834719245500636	0.0	0.0	546.6119690326636;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25264043707101985	0.0	0.0	546.6119690326636;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.1894803278032649	0.0	0.0	546.6119690326636;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29277814738748803	0.0	0.0	546.6119690326636;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.21958361054061604	0.0	0.0	546.6119690326636;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2680156701618669	0.0	0.0	546.6119690326636;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.20101175262140017	0.0	0.0	546.6119690326636;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2785506997132444	0.0	0.0	546.6119690326636;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.20891302478493332	0.0	0.0	546.6119690326636;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28239414142251534	0.0	0.0	546.6119690326636;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.21179560606688652	0.0	0.0	546.6119690326636;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25073261934525626	0.0	0.0	546.6119690326636;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.29185010744510437	0.0	0.0	546.6119690326636;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.2188875805838283	0.0	0.0	546.6119690326636;	%DE1 76 PHS (pump mode)
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

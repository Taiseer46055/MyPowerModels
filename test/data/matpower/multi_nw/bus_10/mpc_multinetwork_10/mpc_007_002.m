function mpc = mpc_007_002
%MPC_007_002	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 7 	Weight: 33
%	Time step: 2

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 33;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	949.2957177457894	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	7609.134967658925	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4103.0379251383765	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2396.483078008019	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	5681.291999883989	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	7493.603715490725	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8824.412237512108	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4904.436150843603	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4085.372161290104	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	7356.015322730983	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.6890022840487227	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.428653445336072	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.5326229435685654	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.4014454123918294	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5019842830118395	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.585944335371216	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.310532244135619	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.4442890779192648	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3055625404097446	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2819859286456645	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.7143567273081577	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.18300541605734	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.1513208503229757	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	275218.6726924449	110087.46907697796	3	0.0	1404.1769014920656	0.0	0.0	373.0496651454373;	%0 CCGT
	2	371274.7263757003	148509.8905502801	3	0.0	1894.258808039287	0.0	0.0	169.88008789943405;	%0 OCGT
	2	2512.7491618410354	1005.0996647364142	3	0.0	491.8279823529135	0.0	0.0	1047.573191409757;	%0 biomass
	2	0.0	0.0	3	0.0	0.8153730553782292	0.0	0.0	724.8496179704457;	%0 offwind
	2	0.0	0.0	3	0.0	0.8110488807463043	0.0	0.0	361.9673863765854;	%0 onwind
	2	0.0	0.0	3	0.0	0.6490295407330537	0.0	0.0	134.11739167186607;	%0 solar
	2	276454.9892289165	110581.99569156661	3	0.0	1410.4846389230436	0.0	0.0	373.0496651454373;	%1 CCGT
	2	355576.7039806721	142230.68159226884	3	0.0	1814.1668570442453	0.0	0.0	169.88008789943405;	%1 OCGT
	2	2512.69624833307	1005.0784993332279	3	0.0	491.81762543219213	0.0	0.0	1047.573191409757;	%1 biomass
	2	0.0	0.0	3	0.0	0.825481003344585	0.0	0.0	361.9673863765854;	%1 onwind
	2	0.0	0.0	3	0.0	0.6643114613186722	0.0	0.0	134.11739167186607;	%1 solar
	2	263545.2946773218	105418.1178709287	3	0.0	1344.6188503944986	0.0	0.0	373.0496651454373;	%2 CCGT
	2	365350.92647683085	146140.37059073235	3	0.0	1864.0353391675042	0.0	0.0	169.88008789943405;	%2 OCGT
	2	2512.6755288876398	1005.0702115550561	3	0.0	491.81356995256215	0.0	0.0	1047.5731914097569;	%2 biomass
	2	0.0	0.0	3	0.0	0.8287295061754857	0.0	0.0	361.9673863765854;	%2 onwind
	2	0.0	0.0	3	0.0	0.6587063790987753	0.0	0.0	134.11739167186607;	%2 solar
	2	244031.33846355663	97612.53538542266	3	0.0	1245.0578493038604	0.0	0.0	373.0496651454373;	%3 CCGT
	2	350560.8298582549	140224.33194330195	3	0.0	1788.5756625421163	0.0	0.0	169.88008789943405;	%3 OCGT
	2	2512.74706544251	1005.098826177004	3	0.0	491.82757201849864	0.0	0.0	1047.573191409757;	%3 biomass
	2	0.0	0.0	3	0.0	0.8297807104355084	0.0	0.0	742.1837821731442;	%3 offwind
	2	0.0	0.0	3	0.0	0.8149317869380144	0.0	0.0	361.96738637658547;	%3 onwind
	2	0.0	0.0	3	0.0	0.6550604577109422	0.0	0.0	134.11739167186607;	%3 solar
	2	261305.5149856772	104522.20599427087	3	0.0	1333.1914029881489	0.0	0.0	373.0496651454373;	%4 CCGT
	2	360641.7002325025	144256.680093001	3	0.0	1840.0086746556249	0.0	0.0	169.88008789943407;	%4 OCGT
	2	2512.749788653332	1005.0999154613328	3	0.0	491.8281050407774	0.0	0.0	1047.573191409757;	%4 biomass
	2	0.0	0.0	3	0.0	0.8171060117131382	0.0	0.0	361.9673863765854;	%4 onwind
	2	0.0	0.0	3	0.0	0.6453282342924473	0.0	0.0	134.11739167186607;	%4 solar
	2	272066.23541661614	108826.49416664644	3	0.0	1388.093037839878	0.0	0.0	373.0496651454373;	%5 CCGT
	2	362313.79421035305	144925.51768414123	3	0.0	1848.5397663793524	0.0	0.0	169.88008789943405;	%5 OCGT
	2	2512.76188598497	1005.104754393988	3	0.0	491.8304728880348	0.0	0.0	1047.573191409757;	%5 biomass
	2	0.0	0.0	3	0.0	0.8255812770239439	0.0	0.0	361.9673863765854;	%5 onwind
	2	0.0	0.0	3	0.0	0.6679718346363013	0.0	0.0	134.11739167186607;	%5 solar
	2	258772.03226510406	103508.81290604162	3	0.0	1320.2654707403267	0.0	0.0	373.0496651454373;	%6 CCGT
	2	365523.03084054904	146209.2123362196	3	0.0	1864.913422655862	0.0	0.0	169.88008789943407;	%6 OCGT
	2	2512.72669732424	1005.090678929696	3	0.0	491.82358530519474	0.0	0.0	1047.573191409757;	%6 biomass
	2	0.0	0.0	3	0.0	0.8371584938037939	0.0	0.0	361.9673863765854;	%6 onwind
	2	0.0	0.0	3	0.0	0.6666562726708487	0.0	0.0	134.11739167186605;	%6 solar
	2	262983.34605532314	105193.33842212927	3	0.0	1341.7517655883833	0.0	0.0	373.0496651454373;	%7 CCGT
	2	367961.9051669109	147184.76206676435	3	0.0	1877.3566590148514	0.0	0.0	169.88008789943405;	%7 OCGT
	2	2512.7596306940204	1005.1038522776082	3	0.0	491.8300314531259	0.0	0.0	1047.573191409757;	%7 biomass
	2	0.0	0.0	3	0.0	0.8306820034305787	0.0	0.0	361.9673863765854;	%7 onwind
	2	0.0	0.0	3	0.0	0.6549314766627082	0.0	0.0	134.11739167186607;	%7 solar
	2	282835.3919062249	113134.15676248998	3	0.0	1443.03771380727	0.0	0.0	373.0496651454373;	%8 CCGT
	2	369985.7691755594	147994.30767022376	3	0.0	1887.6824957936703	0.0	0.0	169.88008789943405;	%8 OCGT
	2	2512.779216313493	1005.1116865253972	3	0.0	491.833865005577	0.0	0.0	1047.573191409757;	%8 biomass
	2	0.0	0.0	3	0.0	0.8362743924922019	0.0	0.0	361.9673863765854;	%8 onwind
	2	0.0	0.0	3	0.0	0.6583358319864403	0.0	0.0	134.11739167186607;	%8 solar
	2	284835.3872310395	113934.1548924158	3	0.0	1453.2417715869362	0.0	0.0	373.0496651454373;	%9 CCGT
	2	371378.2520161492	148551.3008064597	3	0.0	1894.787000082394	0.0	0.0	169.88008789943407;	%9 OCGT
	2	2512.6803252703335	1005.0721301081334	3	0.0	491.8145087630325	0.0	0.0	1047.573191409757;	%9 biomass
	2	0.0	0.0	3	0.0	0.817552134241991	0.0	0.0	766.5235539149801;	%9 offwind
	2	0.0	0.0	3	0.0	0.8242850440136428	0.0	0.0	361.9673863765854;	%9 onwind
	2	0.0	0.0	3	0.0	0.6597158418544935	0.0	0.0	134.11739167186607;	%9 solar
	2	0.0	0.0	3	0.0	0.31904686857606035	0.0	0.0	5852392.148443051;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.23928515143204526	0.0	0.0	5852392.148443051;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34759548085118835	0.0	0.0	5852392.148443051;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.2606966106383913	0.0	0.0	5852392.148443051;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3626972450121839	0.0	0.0	5852392.148443051;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.27202293375913794	0.0	0.0	5852392.148443051;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3028886551525132	0.0	0.0	5852392.148443051;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.22716649136438488	0.0	0.0	5852392.148443051;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3603496610707074	0.0	0.0	5852392.148443051;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.2702622458030306	0.0	0.0	5852392.148443051;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34665661785790935	0.0	0.0	5852392.148443051;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.259992463393432	0.0	0.0	5852392.148443051;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32997108482855814	0.0	0.0	5852392.148443051;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.2474783136214186	0.0	0.0	5852392.148443051;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33494537941096664	0.0	0.0	5852392.148443051;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.33117545729982906	0.0	0.0	5852392.148443051;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.24838159297487178	0.0	0.0	5852392.148443051;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32249036567090483	0.0	0.0	5852392.148443051;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.24186777425317862	0.0	0.0	5852392.148443051;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31656288905563723	0.0	0.0	5852392.148443051;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.23742216679172792	0.0	0.0	5852392.148443051;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3031553345845837	0.0	0.0	5852392.148443051;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.2273665009384378	0.0	0.0	5852392.148443051;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3343709187781376	0.0	0.0	5852392.148443051;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.2507781890836032	0.0	0.0	5852392.148443051;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3069361654822326	0.0	0.0	5852392.148443051;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.23020212411167446	0.0	0.0	5852392.148443051;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30878275642013536	0.0	0.0	5852392.148443051;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.23158706731510154	0.0	0.0	5852392.148443051;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.35783995791804096	0.0	0.0	5852392.148443051;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.26837996843853074	0.0	0.0	5852392.148443051;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32757470797561505	0.0	0.0	5852392.148443051;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.2456810309817113	0.0	0.0	5852392.148443051;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34045085520507656	0.0	0.0	5852392.148443051;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.2553381414038074	0.0	0.0	5852392.148443051;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3451483950719632	0.0	0.0	5852392.148443051;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.2588612963039724	0.0	0.0	5852392.148443051;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3064509791997576	0.0	0.0	5852392.148443051;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.35670568687734977	0.0	0.0	5852392.148443051;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.26752926515801234	0.0	0.0	5852392.148443051;	%DE1 76 PHS (pump mode)
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

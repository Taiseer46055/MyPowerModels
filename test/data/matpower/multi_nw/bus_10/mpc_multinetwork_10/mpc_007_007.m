function mpc = mpc_007_007
%MPC_007_007	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 7 	Weight: 28
%	Time step: 7

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 28;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	851.949130195887	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	9347.04189416738	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6930.70105733372	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4606.682034853788	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2380.6012772586005	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	2000.8156120965327	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	12378.10515040031	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5611.741865431444	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3875.8337070251787	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	9004.723624816465	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.008648265253569867	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.64953738446101	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.04021691779634977	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.867544322828569	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.4600940388055825	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.9039472285973598	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.5754230814624142	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.5777294304463245	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.1782835200001704	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.619721569763717	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.2404299119954105	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.707820590077878	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.10447015594293832	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.475367536497181	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.07163784444302787	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.379855920666659	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.328095084806181	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6788996592679277	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.663022074194137	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	233518.8737996502	93407.54951986008	3	0.0	1191.4228255084192	0.0	0.0	70902.04544824797;	%0 CCGT
	2	315020.9799551397	126008.39198205585	3	0.0	1607.2498977303042	0.0	0.0	32287.512463795472;	%0 OCGT
	2	2132.029591865121	852.8118367460484	3	0.0	417.30859108732056	0.0	0.0	8888.499805900969;	%0 biomass
	2	0.0	0.0	3	0.0	0.6918316833512248	0.0	0.0	10455.406610725218;	%0 offwind
	2	0.0	0.0	3	0.0	0.688162686693834	0.0	0.0	11670.706033475357;	%0 onwind
	2	0.0	0.0	3	0.0	0.550691731531076	0.0	0.0	3186.304093052818;	%0 solar
	2	234567.86964877768	93827.14785951107	3	0.0	1196.7748451468249	0.0	0.0	70902.04544824797;	%1 CCGT
	2	301701.4458017824	120680.57832071296	3	0.0	1539.2930908254202	0.0	0.0	32287.512463795472;	%1 OCGT
	2	2131.984695555332	852.7938782221328	3	0.0	417.2998033970115	0.0	0.0	8888.499805900969;	%1 biomass
	2	0.0	0.0	3	0.0	0.7004081240499509	0.0	0.0	11670.706033475357;	%1 onwind
	2	0.0	0.0	3	0.0	0.5636582096037219	0.0	0.0	3186.304093052818;	%1 solar
	2	223614.18942318214	89445.67576927284	3	0.0	1140.8887215468474	0.0	0.0	70902.04544824797;	%2 CCGT
	2	309994.7254954928	123997.89019819714	3	0.0	1581.6057423239429	0.0	0.0	32287.512463795472;	%2 OCGT
	2	2131.967115419816	852.7868461679263	3	0.0	417.2963623839921	0.0	0.0	8888.499805900969;	%2 biomass
	2	0.0	0.0	3	0.0	0.7031644294822303	0.0	0.0	11670.706033475357;	%2 onwind
	2	0.0	0.0	3	0.0	0.5589023822656275	0.0	0.0	3186.304093052818;	%2 solar
	2	207056.89324180564	82822.75729672225	3	0.0	1056.4127206214573	0.0	0.0	70902.04544824797;	%3 CCGT
	2	297445.5526070041	118978.22104280164	3	0.0	1517.579350035735	0.0	0.0	32287.512463795472;	%3 OCGT
	2	2132.0278131027358	852.8111252410943	3	0.0	417.3082429247867	0.0	0.0	8888.499805900969;	%3 biomass
	2	0.0	0.0	3	0.0	0.7040563603695222	0.0	0.0	10705.438797406565;	%3 offwind
	2	0.0	0.0	3	0.0	0.691457273765588	0.0	0.0	11670.706033475359;	%3 onwind
	2	0.0	0.0	3	0.0	0.5558088732092843	0.0	0.0	3186.304093052818;	%3 solar
	2	221713.77029087764	88685.50811635105	3	0.0	1131.192705565702	0.0	0.0	70902.04544824797;	%4 CCGT
	2	305999.01837909303	122399.6073516372	3	0.0	1561.2194815259847	0.0	0.0	32287.512463795472;	%4 OCGT
	2	2132.0301237058575	852.8120494823429	3	0.0	417.30869518611416	0.0	0.0	8888.499805900969;	%4 biomass
	2	0.0	0.0	3	0.0	0.6933020705444809	0.0	0.0	11670.706033475357;	%4 onwind
	2	0.0	0.0	3	0.0	0.547551229096622	0.0	0.0	3186.304093052818;	%4 solar
	2	230844.07853531063	92337.63141412426	3	0.0	1177.775910894442	0.0	0.0	70902.04544824797;	%5 CCGT
	2	307417.76478454197	122967.1059138168	3	0.0	1568.457983594602	0.0	0.0	32287.512463795472;	%5 OCGT
	2	2132.0403881084594	852.8161552433837	3	0.0	417.3107042686356	0.0	0.0	8888.499805900969;	%5 biomass
	2	0.0	0.0	3	0.0	0.7004932047475887	0.0	0.0	11670.706033475357;	%5 onwind
	2	0.0	0.0	3	0.0	0.5667639809035283	0.0	0.0	3186.304093052818;	%5 solar
	2	219564.14858857315	87825.65943542926	3	0.0	1120.2252479008832	0.0	0.0	70902.04544824797;	%6 CCGT
	2	310140.75344046584	124056.30137618633	3	0.0	1582.3507828595193	0.0	0.0	32287.512463795472;	%6 OCGT
	2	2132.0105310629915	852.8042124251965	3	0.0	417.3048602589531	0.0	0.0	8888.499805900969;	%6 biomass
	2	0.0	0.0	3	0.0	0.710316297772916	0.0	0.0	11670.706033475357;	%6 onwind
	2	0.0	0.0	3	0.0	0.565647746508599	0.0	0.0	3186.3040930528177;	%6 solar
	2	223137.38453178934	89254.95381271574	3	0.0	1138.4560435295373	0.0	0.0	70902.04544824797;	%7 CCGT
	2	312210.1013537426	124884.04054149703	3	0.0	1592.9086803762375	0.0	0.0	32287.512463795472;	%7 OCGT
	2	2132.03847452826	852.8153898113039	3	0.0	417.3103297178038	0.0	0.0	8888.499805900969;	%7 biomass
	2	0.0	0.0	3	0.0	0.704821093819885	0.0	0.0	11670.706033475357;	%7 onwind
	2	0.0	0.0	3	0.0	0.5556994347441161	0.0	0.0	3186.304093052818;	%7 solar
	2	239981.54464770597	95992.61785908241	3	0.0	1224.3956359576837	0.0	0.0	70902.04544824797;	%8 CCGT
	2	313927.31930047466	125570.92772018987	3	0.0	1601.669996430993	0.0	0.0	32287.512463795472;	%8 OCGT
	2	2132.0550926296305	852.8220370518522	3	0.0	417.3135824289744	0.0	0.0	8888.499805900969;	%8 biomass
	2	0.0	0.0	3	0.0	0.7095661512055047	0.0	0.0	11670.706033475357;	%8 onwind
	2	0.0	0.0	3	0.0	0.5585879786551614	0.0	0.0	3186.304093052818;	%8 solar
	2	241678.51037785172	96671.40415114068	3	0.0	1233.0536243767942	0.0	0.0	70902.04544824797;	%9 CCGT
	2	315108.8198924902	126043.52795699611	3	0.0	1607.6980606759705	0.0	0.0	32287.512463795472;	%9 OCGT
	2	2131.9711850778585	852.7884740311434	3	0.0	417.29715895045183	0.0	0.0	8888.499805900969;	%9 biomass
	2	0.0	0.0	3	0.0	0.6936805987507803	0.0	0.0	11056.52156556153;	%9 offwind
	2	0.0	0.0	3	0.0	0.6993933706782424	0.0	0.0	11670.706033475357;	%9 onwind
	2	0.0	0.0	3	0.0	0.5597588961189642	0.0	0.0	3186.304093052818;	%9 solar
	2	0.0	0.0	3	0.0	0.27070643394332394	0.0	0.0	566.8568567746141;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.20302982545749296	0.0	0.0	566.8568567746141;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2949294989040386	0.0	0.0	566.8568567746141;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.22119712417802895	0.0	0.0	566.8568567746141;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3077431169800348	0.0	0.0	566.8568567746141;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.2308073377350261	0.0	0.0	566.8568567746141;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25699643467485966	0.0	0.0	566.8568567746141;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.19274732600614475	0.0	0.0	566.8568567746141;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30575122757514567	0.0	0.0	566.8568567746141;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.22931342068135924	0.0	0.0	566.8568567746141;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29413288787943825	0.0	0.0	566.8568567746141;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.2205996659095787	0.0	0.0	566.8568567746141;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2799754659151402	0.0	0.0	566.8568567746141;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.20998159943635517	0.0	0.0	566.8568567746141;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2841960795002141	0.0	0.0	566.8568567746141;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.28099735770894585	0.0	0.0	566.8568567746141;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.2107480182817094	0.0	0.0	566.8568567746141;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2736281890541011	0.0	0.0	566.8568567746141;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.20522114179057582	0.0	0.0	566.8568567746141;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2685988149562983	0.0	0.0	566.8568567746141;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.2014491112172237	0.0	0.0	566.8568567746141;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.257222708132374	0.0	0.0	566.8568567746141;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.19291703109928052	0.0	0.0	566.8568567746141;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28370865835720765	0.0	0.0	566.8568567746141;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.21278149376790573	0.0	0.0	566.8568567746141;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26043068586371254	0.0	0.0	566.8568567746141;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.1953230143977844	0.0	0.0	566.8568567746141;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26199749029587244	0.0	0.0	566.8568567746141;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.19649811772190434	0.0	0.0	566.8568567746141;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30362178247591354	0.0	0.0	566.8568567746141;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.22771633685693515	0.0	0.0	566.8568567746141;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.27794217646415825	0.0	0.0	566.8568567746141;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.20845663234811868	0.0	0.0	566.8568567746141;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28886739229521646	0.0	0.0	566.8568567746141;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.21665054422141233	0.0	0.0	566.8568567746141;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29285318369742336	0.0	0.0	566.8568567746141;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.21963988777306753	0.0	0.0	566.8568567746141;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2600190126543398	0.0	0.0	566.8568567746141;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.30265937068381193	0.0	0.0	566.8568567746141;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.22699452801285896	0.0	0.0	566.8568567746141;	%DE1 76 PHS (pump mode)
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

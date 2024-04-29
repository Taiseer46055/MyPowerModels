function mpc = mpc_008_002
%MPC_008_002	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 8 	Weight: 42
%	Time step: 2

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 42;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	-7.7114259794085696	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4688.12367309191	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	2897.2226197570603	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	754.4272135726462	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	1253.2567595639218	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-612.586346156596	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8189.822545179656	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3439.625848683609	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	1487.4377684929116	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	5944.928135851088	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.9368862200011119	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.131621550354462	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.1607279416304257	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.201706298185644	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.383268278722735	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5823421684504106	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.152241303525574	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3038460222496533	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2833773530901476	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.905373084457824	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.071227084772786	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.4984702435669086	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.578753957691311	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	350278.3106994753	140111.3242797901	3	0.0	1787.134238262629	0.0	0.0	106353.06817237195;	%0 CCGT
	2	472531.46993270953	189012.58797308378	3	0.0	2410.874846595456	0.0	0.0	48431.26869569321;	%0 OCGT
	2	3198.0443877976813	1279.2177551190725	3	0.0	625.9628866309808	0.0	0.0	13332.749708851454;	%0 biomass
	2	0.0	0.0	3	0.0	1.037747525026837	0.0	0.0	15683.109916087826;	%0 offwind
	2	0.0	0.0	3	0.0	1.032244030040751	0.0	0.0	17506.059050213036;	%0 onwind
	2	0.0	0.0	3	0.0	0.8260375972966139	0.0	0.0	4779.456139579227;	%0 solar
	2	351851.8044731665	140740.7217892666	3	0.0	1795.1622677202372	0.0	0.0	106353.06817237195;	%1 CCGT
	2	452552.16870267363	181020.86748106944	3	0.0	2308.9396362381303	0.0	0.0	48431.26869569321;	%1 OCGT
	2	3197.9770433329977	1279.190817333199	3	0.0	625.9497050955173	0.0	0.0	13332.749708851454;	%1 biomass
	2	0.0	0.0	3	0.0	1.0506121860749262	0.0	0.0	17506.059050213036;	%1 onwind
	2	0.0	0.0	3	0.0	0.8454873144055828	0.0	0.0	4779.456139579227;	%1 solar
	2	335421.2841347732	134168.51365390926	3	0.0	1711.333082320271	0.0	0.0	106353.06817237195;	%2 CCGT
	2	464992.08824323927	185996.8352972957	3	0.0	2372.4086134859144	0.0	0.0	48431.26869569321;	%2 OCGT
	2	3197.9506731297233	1279.1802692518895	3	0.0	625.9445435759882	0.0	0.0	13332.749708851454;	%2 biomass
	2	0.0	0.0	3	0.0	1.0547466442233455	0.0	0.0	17506.059050213036;	%2 onwind
	2	0.0	0.0	3	0.0	0.8383535733984413	0.0	0.0	4779.456139579227;	%2 solar
	2	310585.33986270847	124234.13594508338	3	0.0	1584.6190809321859	0.0	0.0	106353.06817237195;	%3 CCGT
	2	446168.3289105062	178467.33156420247	3	0.0	2276.3690250536024	0.0	0.0	48431.26869569321;	%3 OCGT
	2	3198.041719654104	1279.2166878616413	3	0.0	625.9623643871802	0.0	0.0	13332.749708851454;	%3 biomass
	2	0.0	0.0	3	0.0	1.0560845405542834	0.0	0.0	16058.158196109847;	%3 offwind
	2	0.0	0.0	3	0.0	1.037185910648382	0.0	0.0	17506.059050213036;	%3 onwind
	2	0.0	0.0	3	0.0	0.8337133098139266	0.0	0.0	4779.456139579227;	%3 solar
	2	332570.65543631645	133028.26217452658	3	0.0	1696.789058348553	0.0	0.0	106353.06817237195;	%4 CCGT
	2	458998.52756863955	183599.4110274558	3	0.0	2341.829222288977	0.0	0.0	48431.26869569321;	%4 OCGT
	2	3198.045185558786	1279.2180742235143	3	0.0	625.9630427791712	0.0	0.0	13332.749708851454;	%4 biomass
	2	0.0	0.0	3	0.0	1.0399531058167213	0.0	0.0	17506.059050213036;	%4 onwind
	2	0.0	0.0	3	0.0	0.8213268436449329	0.0	0.0	4779.456139579227;	%4 solar
	2	346266.117802966	138506.4471211864	3	0.0	1766.663866341663	0.0	0.0	106353.06817237195;	%5 CCGT
	2	461126.647176813	184450.6588707252	3	0.0	2352.686975391903	0.0	0.0	48431.26869569321;	%5 OCGT
	2	3198.0605821626887	1279.2242328650757	3	0.0	625.9660564029534	0.0	0.0	13332.749708851454;	%5 biomass
	2	0.0	0.0	3	0.0	1.0507398071213832	0.0	0.0	17506.059050213036;	%5 onwind
	2	0.0	0.0	3	0.0	0.8501459713552925	0.0	0.0	4779.456139579227;	%5 solar
	2	329346.2228828597	131738.4891531439	3	0.0	1680.337871851325	0.0	0.0	106353.06817237195;	%6 CCGT
	2	465211.1301606988	186084.4520642795	3	0.0	2373.526174289279	0.0	0.0	48431.26869569321;	%6 OCGT
	2	3198.0157965944873	1279.2063186377948	3	0.0	625.9572903884297	0.0	0.0	13332.749708851454;	%6 biomass
	2	0.0	0.0	3	0.0	1.065474446659374	0.0	0.0	17506.059050213036;	%6 onwind
	2	0.0	0.0	3	0.0	0.8484716197628984	0.0	0.0	4779.456139579227;	%6 solar
	2	334706.076797684	133882.4307190736	3	0.0	1707.684065294306	0.0	0.0	106353.06817237195;	%7 CCGT
	2	468315.15203061386	187326.06081224553	3	0.0	2389.3630205643562	0.0	0.0	48431.26869569321;	%7 OCGT
	2	3198.05771179239	1279.2230847169558	3	0.0	625.9654945767057	0.0	0.0	13332.749708851454;	%7 biomass
	2	0.0	0.0	3	0.0	1.0572316407298274	0.0	0.0	17506.059050213036;	%7 onwind
	2	0.0	0.0	3	0.0	0.8335491521161741	0.0	0.0	4779.456139579227;	%7 solar
	2	359972.31697155896	143988.92678862362	3	0.0	1836.5934539365255	0.0	0.0	106353.06817237195;	%8 CCGT
	2	470890.978950712	188356.3915802848	3	0.0	2402.504994646489	0.0	0.0	48431.26869569321;	%8 OCGT
	2	3198.082638944446	1279.2330555777783	3	0.0	625.9703736434617	0.0	0.0	13332.749708851454;	%8 biomass
	2	0.0	0.0	3	0.0	1.064349226808257	0.0	0.0	17506.059050213036;	%8 onwind
	2	0.0	0.0	3	0.0	0.8378819679827422	0.0	0.0	4779.456139579227;	%8 solar
	2	362517.76556677755	145007.106226711	3	0.0	1849.5804365651916	0.0	0.0	106353.06817237195;	%9 CCGT
	2	472663.22983873537	189065.29193549417	3	0.0	2411.5470910139557	0.0	0.0	48431.26869569321;	%9 OCGT
	2	3197.956777616788	1279.1827110467152	3	0.0	625.9457384256777	0.0	0.0	13332.749708851454;	%9 biomass
	2	0.0	0.0	3	0.0	1.0405208981261704	0.0	0.0	16584.782348342294;	%9 offwind
	2	0.0	0.0	3	0.0	1.0490900560173635	0.0	0.0	17506.059050213036;	%9 onwind
	2	0.0	0.0	3	0.0	0.8396383441784463	0.0	0.0	4779.456139579227;	%9 solar
	2	0.0	0.0	3	0.0	0.40605965091498586	0.0	0.0	850.2852851619211;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.3045447381862394	0.0	0.0	850.2852851619211;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4423942483560579	0.0	0.0	850.2852851619211;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.3317956862670434	0.0	0.0	850.2852851619211;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46161467547005225	0.0	0.0	850.2852851619211;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.3462110066025392	0.0	0.0	850.2852851619211;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3854946520122895	0.0	0.0	850.2852851619211;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.2891209890092171	0.0	0.0	850.2852851619211;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.45862684136271853	0.0	0.0	850.2852851619211;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.3439701310220389	0.0	0.0	850.2852851619211;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4411993318191573	0.0	0.0	850.2852851619211;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.330899498864368	0.0	0.0	850.2852851619211;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41996319887271033	0.0	0.0	850.2852851619211;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.31497239915453273	0.0	0.0	850.2852851619211;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4262941192503212	0.0	0.0	850.2852851619211;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.4214960365634188	0.0	0.0	850.2852851619211;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.31612202742256407	0.0	0.0	850.2852851619211;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4104422835811516	0.0	0.0	850.2852851619211;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.3078317126858637	0.0	0.0	850.2852851619211;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4028982224344474	0.0	0.0	850.2852851619211;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.3021736668258356	0.0	0.0	850.2852851619211;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3858340621985611	0.0	0.0	850.2852851619211;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.28937554664892084	0.0	0.0	850.2852851619211;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4255629875358115	0.0	0.0	850.2852851619211;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.3191722406518586	0.0	0.0	850.2852851619211;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39064602879556876	0.0	0.0	850.2852851619211;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.29298452159667654	0.0	0.0	850.2852851619211;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39299623544380863	0.0	0.0	850.2852851619211;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.2947471765828565	0.0	0.0	850.2852851619211;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4554326737138703	0.0	0.0	850.2852851619211;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.34157450528540273	0.0	0.0	850.2852851619211;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41691326469623735	0.0	0.0	850.2852851619211;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.312684948522178	0.0	0.0	850.2852851619211;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43330108844282467	0.0	0.0	850.2852851619211;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.3249758163321185	0.0	0.0	850.2852851619211;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.439279775546135	0.0	0.0	850.2852851619211;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.32945983165960124	0.0	0.0	850.2852851619211;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3900285189815097	0.0	0.0	850.2852851619211;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.4539890560257179	0.0	0.0	850.2852851619211;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.3404917920192884	0.0	0.0	850.2852851619211;	%DE1 76 PHS (pump mode)
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

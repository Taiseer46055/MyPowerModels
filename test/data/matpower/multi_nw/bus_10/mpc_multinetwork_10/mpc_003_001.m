function mpc = mpc_003_001
%MPC_003_001	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 3 	Weight: 38
%	Time step: 1

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 38;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1938.3416378788008	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	7429.2203482993755	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4746.021963284383	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3317.7504304710255	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	8122.572372707419	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	10309.140862881439	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	1099.7378853916205	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4092.398867899357	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5650.67441153869	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3715.5193302508583	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.309754304732133	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	21.263957879245382	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.45826222946512	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666665	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.256393947814974	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11000000000001	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.786058476940926	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.793552156938613	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.074949619830317	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.325897205018311	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.36270237147436	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66000000000003	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.393010383496133	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.600596804003118	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.123614313444723	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.77339721896391	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	316918.4715852396	126767.38863409583	3	0.0	1616.9309774757119	0.0	0.0	429.57234168262477;	%0 CCGT
	2	427528.47279626096	171011.3891185044	3	0.0	2181.2677183482697	0.0	0.0	195.61949515692407;	%0 OCGT
	2	2893.4687318169495	1157.38749272678	3	0.0	566.3473736185065	0.0	0.0	1206.2964022294173;	%0 biomass
	2	0.0	0.0	3	0.0	0.9389144274052337	0.0	0.0	834.6753176629375;	%0 offwind
	2	0.0	0.0	3	0.0	0.9339350747987746	0.0	0.0	416.8109297669771;	%0 onwind
	2	0.0	0.0	3	0.0	0.7473673499350316	0.0	0.0	154.43820859184578;	%0 solar
	2	318342.1088090554	127336.84352362217	3	0.0	1624.1944326992623	0.0	0.0	429.57234168262477;	%1 CCGT
	2	409451.96215956187	163780.78486382472	3	0.0	2089.0406232630703	0.0	0.0	195.61949515692407;	%1 OCGT
	2	2893.4078011108077	1157.363120444323	3	0.0	566.3354474673728	0.0	0.0	1206.2964022294173;	%1 biomass
	2	0.0	0.0	3	0.0	0.950553882639219	0.0	0.0	416.8109297669771;	%1 onwind
	2	0.0	0.0	3	0.0	0.7649647130336226	0.0	0.0	154.43820859184578;	%1 solar
	2	303476.39993146143	121390.55997258458	3	0.0	1548.34897924215	0.0	0.0	429.57234168262477;	%2 CCGT
	2	420707.12745816883	168282.85098326753	3	0.0	2146.464936011065	0.0	0.0	195.61949515692407;	%2 OCGT
	2	2893.383942355464	1157.3535769421858	3	0.0	566.3307775211322	0.0	0.0	1206.296402229417;	%2 biomass
	2	0.0	0.0	3	0.0	0.9542945828687411	0.0	0.0	416.8109297669771;	%2 onwind
	2	0.0	0.0	3	0.0	0.758510375931923	0.0	0.0	154.43820859184578;	%2 solar
	2	281005.78368530766	112402.31347412306	3	0.0	1433.7029779862635	0.0	0.0	429.57234168262477;	%3 CCGT
	2	403676.10710950557	161470.44284380224	3	0.0	2059.5719750484977	0.0	0.0	195.61949515692407;	%3 OCGT
	2	2893.466317782284	1157.3865271129137	3	0.0	566.3469011122106	0.0	0.0	1206.2964022294173;	%3 biomass
	2	0.0	0.0	3	0.0	0.9555050605014945	0.0	0.0	854.6358703811964;	%3 offwind
	2	0.0	0.0	3	0.0	0.9384063001104408	0.0	0.0	416.81092976697715;	%3 onwind
	2	0.0	0.0	3	0.0	0.7543120422126002	0.0	0.0	154.43820859184578;	%3 solar
	2	300897.2596804768	120358.90387219071	3	0.0	1535.1901004105957	0.0	0.0	429.57234168262477;	%4 CCGT
	2	415284.38208591193	166113.75283436477	3	0.0	2118.7978677852648	0.0	0.0	195.6194951569241;	%4 OCGT
	2	2893.469453600807	1157.3877814403227	3	0.0	566.3475148954407	0.0	0.0	1206.2964022294173;	%4 biomass
	2	0.0	0.0	3	0.0	0.9409099528817955	0.0	0.0	416.8109297669771;	%4 onwind
	2	0.0	0.0	3	0.0	0.7431052394882727	0.0	0.0	154.43820859184578;	%4 solar
	2	313288.3922979216	125315.35691916864	3	0.0	1598.4101647853142	0.0	0.0	429.57234168262477;	%5 CCGT
	2	417209.8236361641	166883.92945446566	3	0.0	2128.621549164103	0.0	0.0	195.61949515692407;	%5 OCGT
	2	2893.4833838614804	1157.3933535445922	3	0.0	566.350241507434	0.0	0.0	1206.2964022294173;	%5 biomass
	2	0.0	0.0	3	0.0	0.950669349300299	0.0	0.0	416.8109297669771;	%5 onwind
	2	0.0	0.0	3	0.0	0.7691796883690741	0.0	0.0	154.43820859184578;	%5 solar
	2	297979.915941635	119191.96637665399	3	0.0	1520.3056935797701	0.0	0.0	429.57234168262477;	%6 CCGT
	2	420905.3082406322	168362.1232962529	3	0.0	2147.4760624522046	0.0	0.0	195.6194951569241;	%6 OCGT
	2	2893.442863585489	1157.3771454341954	3	0.0	566.3423103514364	0.0	0.0	1206.2964022294173;	%6 biomass
	2	0.0	0.0	3	0.0	0.9640006898346718	0.0	0.0	416.8109297669771;	%6 onwind
	2	0.0	0.0	3	0.0	0.7676647988330986	0.0	0.0	154.43820859184575;	%6 solar
	2	302829.30757885694	121131.72303154279	3	0.0	1545.0474876472292	0.0	0.0	429.57234168262477;	%7 CCGT
	2	423713.7089800792	169485.4835920317	3	0.0	2161.8046376534653	0.0	0.0	195.61949515692407;	%7 OCGT
	2	2893.4807868597813	1157.3923147439125	3	0.0	566.349733188448	0.0	0.0	1206.2964022294173;	%7 biomass
	2	0.0	0.0	3	0.0	0.9565429130412725	0.0	0.0	416.8109297669771;	%7 onwind
	2	0.0	0.0	3	0.0	0.7541635185813004	0.0	0.0	154.43820859184578;	%7 solar
	2	325689.23916474386	130275.69566589755	3	0.0	1661.6797916568564	0.0	0.0	429.57234168262477;	%8 CCGT
	2	426044.2190506442	170417.68762025767	3	0.0	2173.6949951563474	0.0	0.0	195.61949515692407;	%8 OCGT
	2	2893.503339997356	1157.4013359989424	3	0.0	566.3541475821796	0.0	0.0	1206.2964022294173;	%8 biomass
	2	0.0	0.0	3	0.0	0.9629826337788993	0.0	0.0	416.8109297669771;	%8 onwind
	2	0.0	0.0	3	0.0	0.7580836853177191	0.0	0.0	154.43820859184578;	%8 solar
	2	327992.2640842273	131196.90563369094	3	0.0	1673.429918797078	0.0	0.0	429.57234168262477;	%9 CCGT
	2	427647.6841398082	171059.07365592328	3	0.0	2181.875939488817	0.0	0.0	195.6194951569241;	%9 OCGT
	2	2893.389465462808	1157.3557861851232	3	0.0	566.3318585756132	0.0	0.0	1206.2964022294173;	%9 biomass
	2	0.0	0.0	3	0.0	0.9414236697332018	0.0	0.0	882.6634863263406;	%9 offwind
	2	0.0	0.0	3	0.0	0.9491767173490432	0.0	0.0	416.8109297669771;	%9 onwind
	2	0.0	0.0	3	0.0	0.7596727875900229	0.0	0.0	154.43820859184578;	%9 solar
	2	0.0	0.0	3	0.0	0.36738730320879676	0.0	0.0	6739118.231540483;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.2755404774065976	0.0	0.0	6739118.231540483;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4002614627983381	0.0	0.0	6739118.231540483;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.30019609709875356	0.0	0.0	6739118.231540483;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.417651373044333	0.0	0.0	6739118.231540483;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.31323852978324973	0.0	0.0	6739118.231540483;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3487808756301667	0.0	0.0	6739118.231540483;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.261585656722625	0.0	0.0	6739118.231540483;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41494809456626913	0.0	0.0	6739118.231540483;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.31121107092470185	0.0	0.0	6739118.231540483;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39918034783638046	0.0	0.0	6739118.231540483;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.29938526087728534	0.0	0.0	6739118.231540483;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37996670374197605	0.0	0.0	6739118.231540483;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.28497502780648204	0.0	0.0	6739118.231540483;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38569467932171914	0.0	0.0	6739118.231540483;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.38135355689071226	0.0	0.0	6739118.231540483;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.2860151676680342	0.0	0.0	6739118.231540483;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3713525422877086	0.0	0.0	6739118.231540483;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.2785144067157814	0.0	0.0	6739118.231540483;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.36452696315497624	0.0	0.0	6739118.231540483;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.27339522236623215	0.0	0.0	6739118.231540483;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34908796103679335	0.0	0.0	6739118.231540483;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.261815970777595	0.0	0.0	6739118.231540483;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38503317919906754	0.0	0.0	6739118.231540483;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.28877488439930066	0.0	0.0	6739118.231540483;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3534416451007527	0.0	0.0	6739118.231540483;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.26508123382556453	0.0	0.0	6739118.231540483;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3555680225443983	0.0	0.0	6739118.231540483;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.26667601690829873	0.0	0.0	6739118.231540483;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41205813336016833	0.0	0.0	6739118.231540483;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.30904360002012626	0.0	0.0	6739118.231540483;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3772072394870719	0.0	0.0	6739118.231540483;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.28290542961530396	0.0	0.0	6739118.231540483;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3920343181149366	0.0	0.0	6739118.231540483;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.29402573858620246	0.0	0.0	6739118.231540483;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3974436064465031	0.0	0.0	6739118.231540483;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.2980827048348773	0.0	0.0	6739118.231540483;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3528829457451755	0.0	0.0	6739118.231540483;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.4107520030708876	0.0	0.0	6739118.231540483;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.3080640023031657	0.0	0.0	6739118.231540483;	%DE1 76 PHS (pump mode)
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

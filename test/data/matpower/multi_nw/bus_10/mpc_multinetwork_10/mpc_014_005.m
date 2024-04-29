function mpc = mpc_014_005
%MPC_014_005	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 14 	Weight: 23
%	Time step: 5

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 23;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	903.303679963819	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6531.699158838783	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5759.1364611492445	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	1849.9288475017186	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2902.0414761641214	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	1239.5606037602327	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	7258.327522573123	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4067.160414502347	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6564.182846636529	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	6583.915502522449	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.02834928686514	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.421029702705386	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.0929494734058396	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.9469916536248288	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.19298641960951277	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0960318676333243	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8776045582030063	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.203966721618839	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.236674950296527	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.7162114885538289	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.122146192165373	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5312718589913565	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.22994663303354013	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.17216059126433425	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.458750085062257	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.008739994659694881	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.047037844175225	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.993175361090414	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.116107589138666	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.2327882282053373	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.589392312092052	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.056879978412296	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.48314580453822703	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	191819.07490685553	76727.62996274221	3	0.0	978.668749524773	0.0	0.0	58240.96590391797;	%0 CCGT
	2	258767.23353457902	103506.8934138316	3	0.0	1320.2409874213213	0.0	0.0	26521.88523811771;	%0 OCGT
	2	1751.3100218892064	700.5240087556826	3	0.0	342.7891998217276	0.0	0.0	7301.267697704367;	%0 biomass
	2	0.0	0.0	3	0.0	0.5682903113242204	0.0	0.0	8588.369715952856;	%0 offwind
	2	0.0	0.0	3	0.0	0.5652764926413636	0.0	0.0	9586.651384640472;	%0 onwind
	2	0.0	0.0	3	0.0	0.4523539223290981	0.0	0.0	2617.3212192933865;	%0 solar
	2	192680.7500686388	77072.30002745552	3	0.0	983.065051370606	0.0	0.0	58240.96590391797;	%1 CCGT
	2	247826.1876228927	99130.47504915707	3	0.0	1264.4193246065952	0.0	0.0	26521.88523811771;	%1 OCGT
	2	1751.273142777594	700.5092571110376	3	0.0	342.7819813618309	0.0	0.0	7301.267697704367;	%1 biomass
	2	0.0	0.0	3	0.0	0.5753352447553167	0.0	0.0	9586.651384640472;	%1 onwind
	2	0.0	0.0	3	0.0	0.46300495788877155	0.0	0.0	2617.3212192933865;	%1 solar
	2	183683.08416904247	73473.23366761698	3	0.0	937.1585926991961	0.0	0.0	58240.96590391797;	%2 CCGT
	2	254638.52451415482	101855.40980566194	3	0.0	1299.1761454803816	0.0	0.0	26521.88523811771;	%2 OCGT
	2	1751.2587019519915	700.5034807807966	3	0.0	342.7791548154221	0.0	0.0	7301.267697704367;	%2 biomass
	2	0.0	0.0	3	0.0	0.5775993527889749	0.0	0.0	9586.651384640472;	%2 onwind
	2	0.0	0.0	3	0.0	0.45909838543247977	0.0	0.0	2617.3212192933865;	%2 solar
	2	170082.44802005464	68032.97920802185	3	0.0	867.7675919390541	0.0	0.0	58240.96590391797;	%3 CCGT
	2	244330.27535575337	97732.11014230135	3	0.0	1246.5830375293538	0.0	0.0	26521.88523811771;	%3 OCGT
	2	1751.3085607629614	700.5234243051846	3	0.0	342.7889138310748	0.0	0.0	7301.267697704367;	%3 biomass
	2	0.0	0.0	3	0.0	0.5783320103035362	0.0	0.0	8793.753297869678;	%3 offwind
	2	0.0	0.0	3	0.0	0.5679827605931616	0.0	0.0	9586.651384640474;	%3 onwind
	2	0.0	0.0	3	0.0	0.4565572887076264	0.0	0.0	2617.3212192933865;	%3 solar
	2	182122.02559607805	72848.81023843122	3	0.0	929.1940081432552	0.0	0.0	58240.96590391797;	%4 CCGT
	2	251356.33652568355	100542.53461027342	3	0.0	1282.4302883963446	0.0	0.0	26521.88523811771;	%4 OCGT
	2	1751.310458758383	700.5241835033531	3	0.0	342.7892853314509	0.0	0.0	7301.267697704367;	%4 biomass
	2	0.0	0.0	3	0.0	0.5694981293758236	0.0	0.0	9586.651384640472;	%4 onwind
	2	0.0	0.0	3	0.0	0.44977422390079663	0.0	0.0	2617.3212192933865;	%4 solar
	2	189621.9216540052	75848.76866160207	3	0.0	967.4587839490059	0.0	0.0	58240.96590391797;	%5 CCGT
	2	252521.7353587309	101008.69414349236	3	0.0	1288.3762008098518	0.0	0.0	26521.88523811771;	%5 OCGT
	2	1751.3188902319487	700.5275560927795	3	0.0	342.79093564923636	0.0	0.0	7301.267697704367;	%5 biomass
	2	0.0	0.0	3	0.0	0.5754051324712336	0.0	0.0	9586.651384640472;	%5 onwind
	2	0.0	0.0	3	0.0	0.4655561271707554	0.0	0.0	2617.3212192933865;	%5 solar
	2	180356.2649120422	72142.50596481688	3	0.0	920.1850250614399	0.0	0.0	58240.96590391797;	%6 CCGT
	2	254758.47604038264	101903.39041615305	3	0.0	1299.7881430631764	0.0	0.0	26521.88523811771;	%6 OCGT
	2	1751.294364801743	700.5177459206972	3	0.0	342.7861352127115	0.0	0.0	7301.267697704367;	%6 biomass
	2	0.0	0.0	3	0.0	0.5834741017420382	0.0	0.0	9586.651384640472;	%6 onwind
	2	0.0	0.0	3	0.0	0.46463922034634914	0.0	0.0	2617.321219293386;	%6 solar
	2	183291.42300825552	73316.5692033022	3	0.0	935.1603214706914	0.0	0.0	58240.96590391797;	%7 CCGT
	2	256458.29754057425	102583.31901622971	3	0.0	1308.4607017376238	0.0	0.0	26521.88523811771;	%7 OCGT
	2	1751.3173183624992	700.5269273449996	3	0.0	342.7906279824817	0.0	0.0	7301.267697704367;	%7 biomass
	2	0.0	0.0	3	0.0	0.5789601842091913	0.0	0.0	9586.651384640472;	%7 onwind
	2	0.0	0.0	3	0.0	0.4564673928255239	0.0	0.0	2617.3212192933865;	%7 solar
	2	197127.69738918706	78851.07895567484	3	0.0	1005.7535581080973	0.0	0.0	58240.96590391797;	%8 CCGT
	2	257868.86942538992	103147.54777015596	3	0.0	1315.6574970683157	0.0	0.0	26521.88523811771;	%8 OCGT
	2	1751.330968945768	700.5323875783072	3	0.0	342.7932998523719	0.0	0.0	7301.267697704367;	%8 biomass
	2	0.0	0.0	3	0.0	0.5828579099188074	0.0	0.0	9586.651384640472;	%8 onwind
	2	0.0	0.0	3	0.0	0.45884012532388263	0.0	0.0	2617.3212192933865;	%8 solar
	2	198521.6335246639	79408.65340986557	3	0.0	1012.8654771666525	0.0	0.0	58240.96590391797;	%9 CCGT
	2	258839.38776883128	103535.75510753252	3	0.0	1320.6091212695471	0.0	0.0	26521.88523811771;	%9 OCGT
	2	1751.2620448853838	700.5048179541535	3	0.0	342.77980913787115	0.0	0.0	7301.267697704367;	%9 biomass
	2	0.0	0.0	3	0.0	0.5698090632595695	0.0	0.0	9082.1427145684;	%9 offwind
	2	0.0	0.0	3	0.0	0.574501697342842	0.0	0.0	9586.651384640472;	%9 onwind
	2	0.0	0.0	3	0.0	0.45980195038343485	0.0	0.0	2617.3212192933865;	%9 solar
	2	0.0	0.0	3	0.0	0.2223659993105875	0.0	0.0	465.63241806486155;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.16677449948294062	0.0	0.0	465.63241806486155;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24226351695688886	0.0	0.0	465.63241806486155;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.18169763771766664	0.0	0.0	465.63241806486155;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25278898894788576	0.0	0.0	465.63241806486155;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.1895917417109143	0.0	0.0	465.63241806486155;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21110421419720615	0.0	0.0	465.63241806486155;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.15832816064790461	0.0	0.0	465.63241806486155;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25115279407958396	0.0	0.0	465.63241806486155;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.188364595559688	0.0	0.0	465.63241806486155;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24160915790096713	0.0	0.0	465.63241806486155;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.18120686842572534	0.0	0.0	465.63241806486155;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22997984700172233	0.0	0.0	465.63241806486155;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.17248488525129174	0.0	0.0	465.63241806486155;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23344677958946158	0.0	0.0	465.63241806486155;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.23081925811806267	0.0	0.0	465.63241806486155;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.173114443588547	0.0	0.0	465.63241806486155;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2247660124372973	0.0	0.0	465.63241806486155;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.16857450932797297	0.0	0.0	465.63241806486155;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22063474085695928	0.0	0.0	465.63241806486155;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.16547605564271944	0.0	0.0	465.63241806486155;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2112900816801644	0.0	0.0	465.63241806486155;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.1584675612601233	0.0	0.0	465.63241806486155;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23304639793627774	0.0	0.0	465.63241806486155;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.1747847984522083	0.0	0.0	465.63241806486155;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21392520624519243	0.0	0.0	465.63241806486155;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.16044390468389433	0.0	0.0	465.63241806486155;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21521222417160948	0.0	0.0	465.63241806486155;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.16140916812870712	0.0	0.0	465.63241806486155;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2494036070337861	0.0	0.0	465.63241806486155;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.18705270527533957	0.0	0.0	465.63241806486155;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2283096449527014	0.0	0.0	465.63241806486155;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.17123223371452606	0.0	0.0	465.63241806486155;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23728392938535636	0.0	0.0	465.63241806486155;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.17796294703901727	0.0	0.0	465.63241806486155;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24055797232288345	0.0	0.0	465.63241806486155;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.1804184792421626	0.0	0.0	465.63241806486155;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21358704610892199	0.0	0.0	465.63241806486155;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.24861305449027407	0.0	0.0	465.63241806486155;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.18645979086770556	0.0	0.0	465.63241806486155;	%DE1 76 PHS (pump mode)
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

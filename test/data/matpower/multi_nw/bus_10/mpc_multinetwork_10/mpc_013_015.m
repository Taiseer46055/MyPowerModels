function mpc = mpc_013_015
%MPC_013_015	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 13 	Weight: 8
%	Time step: 15

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 8;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	3907.429310821908	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	8977.665191647664	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	7766.9287610681395	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	6167.949021022499	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	8272.476382148083	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	11416.87550255468	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	4489.1380481546985	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4829.2938192074025	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	11456.017820592453	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	7281.966690232101	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.734347044315555	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.93955684794846	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.124156563614549	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.50098956434196	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.698184990965171	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.89429105935477	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.793078812620968	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.77791554193092	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	34.03285011984787	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.483559117781582	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.885540371954532	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.889895219172828	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.96378212113948	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.177997049940551	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.915534852871605	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.045738385090442	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	24.189465984071532	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.955507518659548	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.51502100120281	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.049988837210346	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.027179102330992	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	34.766760911346395	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.735178702198664	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	66719.67822847149	26687.871291388594	3	0.0	340.4065215738341	0.0	0.0	20257.72727092799;	%0 CCGT
	2	90005.99427289705	36002.39770915882	3	0.0	459.2142564943726	0.0	0.0	9225.00356108442;	%0 OCGT
	2	609.1513119614631	243.66052478458525	3	0.0	119.23102602494873	0.0	0.0	2539.5713731145624;	%0 biomass
	2	0.0	0.0	3	0.0	0.19766619524320708	0.0	0.0	2987.2590316357764;	%0 offwind
	2	0.0	0.0	3	0.0	0.19661791048395255	0.0	0.0	3334.4874381358163;	%0 onwind
	2	0.0	0.0	3	0.0	0.15734049472316455	0.0	0.0	910.372598015091;	%0 solar
	2	67019.39132822219	26807.756531288876	3	0.0	341.93567004194995	0.0	0.0	20257.72727092799;	%1 CCGT
	2	86200.41308622355	34480.16523448942	3	0.0	439.79802595012006	0.0	0.0	9225.00356108442;	%1 OCGT
	2	609.1384844443805	243.6553937777522	3	0.0	119.228515256289	0.0	0.0	2539.5713731145624;	%1 biomass
	2	0.0	0.0	3	0.0	0.20011660687141453	0.0	0.0	3334.4874381358163;	%1 onwind
	2	0.0	0.0	3	0.0	0.16104520274392053	0.0	0.0	910.372598015091;	%1 solar
	2	63889.768406623465	25555.907362649385	3	0.0	325.9682061562421	0.0	0.0	20257.72727092799;	%2 CCGT
	2	88569.92157014081	35427.968628056326	3	0.0	451.88735494969796	0.0	0.0	9225.00356108442;	%2 OCGT
	2	609.1334615485188	243.65338461940752	3	0.0	119.22753210971203	0.0	0.0	2539.5713731145624;	%2 biomass
	2	0.0	0.0	3	0.0	0.20090412270920865	0.0	0.0	3334.4874381358163;	%2 onwind
	2	0.0	0.0	3	0.0	0.15968639493303644	0.0	0.0	910.372598015091;	%2 solar
	2	59159.11235480161	23663.644941920644	3	0.0	301.8322058918449	0.0	0.0	20257.72727092799;	%3 CCGT
	2	84984.44360200118	33993.77744080047	3	0.0	433.59410001021	0.0	0.0	9225.00356108442;	%3 OCGT
	2	609.1508037436388	243.6603214974555	3	0.0	119.23092654993907	0.0	0.0	2539.5713731145624;	%3 biomass
	2	0.0	0.0	3	0.0	0.2011589601055778	0.0	0.0	3058.6967992590185;	%3 offwind
	2	0.0	0.0	3	0.0	0.19755922107588228	0.0	0.0	3334.487438135817;	%3 onwind
	2	0.0	0.0	3	0.0	0.15880253520265267	0.0	0.0	910.372598015091;	%3 solar
	2	63346.79151167932	25338.71660467173	3	0.0	323.1979158759149	0.0	0.0	20257.72727092799;	%4 CCGT
	2	87428.29096545515	34971.31638618206	3	0.0	446.0627090074242	0.0	0.0	9225.00356108442;	%4 OCGT
	2	609.1514639159593	243.6605855663837	3	0.0	119.23105576746119	0.0	0.0	2539.5713731145624;	%4 biomass
	2	0.0	0.0	3	0.0	0.19808630586985168	0.0	0.0	3334.4874381358163;	%4 onwind
	2	0.0	0.0	3	0.0	0.15644320831332056	0.0	0.0	910.372598015091;	%4 solar
	2	65955.45101008876	26382.180404035502	3	0.0	336.5074031126977	0.0	0.0	20257.72727092799;	%5 CCGT
	2	87833.6470812977	35133.458832519085	3	0.0	448.1308524556006	0.0	0.0	9225.00356108442;	%5 OCGT
	2	609.1543966024169	243.66175864096678	3	0.0	119.23162979103874	0.0	0.0	2539.5713731145624;	%5 biomass
	2	0.0	0.0	3	0.0	0.2001409156421682	0.0	0.0	3334.4874381358163;	%5 onwind
	2	0.0	0.0	3	0.0	0.16193256597243666	0.0	0.0	910.372598015091;	%5 solar
	2	62732.61388244947	25093.045552979787	3	0.0	320.0643565431095	0.0	0.0	20257.72727092799;	%6 CCGT
	2	88611.6438401331	35444.65753605324	3	0.0	452.10022367414837	0.0	0.0	9225.00356108442;	%6 OCGT
	2	609.1458660179976	243.65834640719902	3	0.0	119.2299600739866	0.0	0.0	2539.5713731145624;	%6 biomass
	2	0.0	0.0	3	0.0	0.20294751364940458	0.0	0.0	3334.4874381358163;	%6 onwind
	2	0.0	0.0	3	0.0	0.1616136418595997	0.0	0.0	910.3725980150908;	%6 solar
	2	63753.5384376541	25501.41537506164	3	0.0	325.2731552941535	0.0	0.0	20257.72727092799;	%7 CCGT
	2	89202.8861010693	35681.15444042772	3	0.0	455.11676582178217	0.0	0.0	9225.00356108442;	%7 OCGT
	2	609.1538498652171	243.66153994608683	3	0.0	119.23152277651538	0.0	0.0	2539.5713731145624;	%7 biomass
	2	0.0	0.0	3	0.0	0.20137745537711	0.0	0.0	3334.4874381358163;	%7 onwind
	2	0.0	0.0	3	0.0	0.15877126706974745	0.0	0.0	910.372598015091;	%7 solar
	2	68566.15561363028	27426.462245452116	3	0.0	349.8273245593382	0.0	0.0	20257.72727092799;	%8 CCGT
	2	89693.51980013562	35877.40792005425	3	0.0	457.6199989802837	0.0	0.0	9225.00356108442;	%8 OCGT
	2	609.1585978941802	243.66343915767206	3	0.0	119.23245212256413	0.0	0.0	2539.5713731145624;	%8 biomass
	2	0.0	0.0	3	0.0	0.20273318605871563	0.0	0.0	3334.4874381358163;	%8 onwind
	2	0.0	0.0	3	0.0	0.15959656533004613	0.0	0.0	910.372598015091;	%8 solar
	2	69051.00296510049	27620.401186040195	3	0.0	352.30103553622695	0.0	0.0	20257.72727092799;	%9 CCGT
	2	90031.09139785435	36012.436559141745	3	0.0	459.3423030502773	0.0	0.0	9225.00356108442;	%9 OCGT
	2	609.1346243079596	243.65384972318384	3	0.0	119.2277597001291	0.0	0.0	2539.5713731145624;	%9 biomass
	2	0.0	0.0	3	0.0	0.19819445678593722	0.0	0.0	3159.0061615890086;	%9 offwind
	2	0.0	0.0	3	0.0	0.19982667733664067	0.0	0.0	3334.4874381358163;	%9 onwind
	2	0.0	0.0	3	0.0	0.1599311131768469	0.0	0.0	910.372598015091;	%9 solar
	2	0.0	0.0	3	0.0	0.07734469541237826	0.0	0.0	161.95910193560402;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.058008521559283696	0.0	0.0	161.95910193560402;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.0842655711154396	0.0	0.0	161.95910193560402;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.0631991783365797	0.0	0.0	161.95910193560402;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08792660485143852	0.0	0.0	161.95910193560402;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.06594495363857888	0.0	0.0	161.95910193560402;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07342755276424562	0.0	0.0	161.95910193560402;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.055070664573184214	0.0	0.0	161.95910193560402;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08735749359289877	0.0	0.0	161.95910193560402;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.06551812019467407	0.0	0.0	161.95910193560402;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08403796796555378	0.0	0.0	161.95910193560402;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.06302847597416533	0.0	0.0	161.95910193560402;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07999299026146864	0.0	0.0	161.95910193560402;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.05999474269610148	0.0	0.0	161.95910193560402;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08119887985720403	0.0	0.0	161.95910193560402;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.0802849593454131	0.0	0.0	161.95910193560402;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.06021371950905983	0.0	0.0	161.95910193560402;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07817948258688602	0.0	0.0	161.95910193560402;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.05863461194016452	0.0	0.0	161.95910193560402;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07674251855894236	0.0	0.0	161.95910193560402;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.05755688891920677	0.0	0.0	161.95910193560402;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07349220232353544	0.0	0.0	161.95910193560402;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.055119151742651584	0.0	0.0	161.95910193560402;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.0810596166734879	0.0	0.0	161.95910193560402;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.06079471250511593	0.0	0.0	161.95910193560402;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07440876738963215	0.0	0.0	161.95910193560402;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.05580657554222411	0.0	0.0	161.95910193560402;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.0748564257988207	0.0	0.0	161.95910193560402;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.056142319349115524	0.0	0.0	161.95910193560402;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08674908070740386	0.0	0.0	161.95910193560402;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.0650618105305529	0.0	0.0	161.95910193560402;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07941205041833092	0.0	0.0	161.95910193560402;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.05955903781374819	0.0	0.0	161.95910193560402;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08253354065577613	0.0	0.0	161.95910193560402;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.0619001554918321	0.0	0.0	161.95910193560402;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.08367233819926381	0.0	0.0	161.95910193560402;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.06275425364944787	0.0	0.0	161.95910193560402;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.07429114647266852	0.0	0.0	161.95910193560402;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.08647410590966055	0.0	0.0	161.95910193560402;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.06485557943224542	0.0	0.0	161.95910193560402;	%DE1 76 PHS (pump mode)
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

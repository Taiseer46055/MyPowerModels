function mpc = mpc_009_008
%MPC_009_008	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 9 	Weight: 17
%	Time step: 8

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 17;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	2296.9827120423583	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	11530.46766512611	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6753.987804616966	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4923.69711548696	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	9625.520451343038	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	14345.186513805696	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3218.7292851965594	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5805.379122034438	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5881.811464535514	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	6222.77384098387	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.228770732459443	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.99362703226669	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0313155735444153	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.37621268427858	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0805413764487803	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.161081411209988	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.4959071895756733	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.702449115572874	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	33.304302183510565	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.9940849631966211	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.567518408502163	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8825992650456027	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.527782089985603	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.812874277599896	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	29.28909660407509	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.3917978094413674	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	21.162044052200706	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2884553939734973	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.967491470866346	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3816652646048757	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.026576335832926	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	29.750798271059754	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.256182542417954	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	141779.3162355019	56711.72649420076	3	0.0	723.3638583443974	0.0	0.0	43047.67045072198;	%0 CCGT
	2	191262.73782990623	76505.09513196249	3	0.0	975.8302950505417	0.0	0.0	19603.132567304394;	%0 OCGT
	2	1294.446537918109	517.7786151672436	3	0.0	253.36593030301606	0.0	0.0	5396.589167868446;	%0 biomass
	2	0.0	0.0	3	0.0	0.4200406648918151	0.0	0.0	6347.925442226025;	%0 offwind
	2	0.0	0.0	3	0.0	0.41781305977839917	0.0	0.0	7085.78580603861;	%0 onwind
	2	0.0	0.0	3	0.0	0.3343485512867247	0.0	0.0	1934.5417707820682;	%0 solar
	2	142416.20657247215	56966.48262898886	3	0.0	726.6132988391437	0.0	0.0	43047.67045072198;	%1 CCGT
	2	183175.87780822505	73270.35112329	3	0.0	934.5708051440051	0.0	0.0	19603.132567304394;	%1 OCGT
	2	1294.4192794443086	517.7677117777234	3	0.0	253.36059491961413	0.0	0.0	5396.589167868446;	%1 biomass
	2	0.0	0.0	3	0.0	0.42524778960175585	0.0	0.0	7085.78580603861;	%1 onwind
	2	0.0	0.0	3	0.0	0.3422210558308311	0.0	0.0	1934.5417707820682;	%1 solar
	2	135765.75786407487	54306.30314562994	3	0.0	692.6824380820144	0.0	0.0	43047.67045072198;	%2 CCGT
	2	188211.08333654923	75284.4333346197	3	0.0	960.2606292681082	0.0	0.0	19603.132567304394;	%2 OCGT
	2	1294.4086057906025	517.763442316241	3	0.0	253.35850573313806	0.0	0.0	5396.589167868446;	%2 biomass
	2	0.0	0.0	3	0.0	0.4269212607570684	0.0	0.0	7085.78580603861;	%2 onwind
	2	0.0	0.0	3	0.0	0.3393335892327024	0.0	0.0	1934.5417707820682;	%2 solar
	2	125713.11375395342	50285.24550158137	3	0.0	641.3934375201704	0.0	0.0	43047.67045072198;	%3 CCGT
	2	180591.94265425252	72236.777061701	3	0.0	921.3874625216963	0.0	0.0	19603.132567304394;	%3 OCGT
	2	1294.4454579552325	517.7781831820929	3	0.0	253.36571891862053	0.0	0.0	5396.589167868446;	%3 biomass
	2	0.0	0.0	3	0.0	0.4274627902243528	0.0	0.0	6499.730698425415;	%3 offwind
	2	0.0	0.0	3	0.0	0.4198133447862498	0.0	0.0	7085.785806038611;	%3 onwind
	2	0.0	0.0	3	0.0	0.3374553873056369	0.0	0.0	1934.5417707820682;	%3 solar
	2	134611.93196231857	53844.77278492742	3	0.0	686.7955712363191	0.0	0.0	43047.67045072198;	%4 CCGT
	2	185785.11830159219	74314.04732063688	3	0.0	947.8832566407764	0.0	0.0	19603.132567304394;	%4 OCGT
	2	1294.4468608214136	517.7787443285654	3	0.0	253.36599350585502	0.0	0.0	5396.589167868446;	%4 biomass
	2	0.0	0.0	3	0.0	0.4209333999734348	0.0	0.0	7085.78580603861;	%4 onwind
	2	0.0	0.0	3	0.0	0.3324418176658062	0.0	0.0	1934.5417707820682;	%4 solar
	2	140155.3333964386	56062.13335857544	3	0.0	715.0782316144827	0.0	0.0	43047.67045072198;	%5 CCGT
	2	186646.50004775764	74658.60001910306	3	0.0	952.2780614681512	0.0	0.0	19603.132567304394;	%5 OCGT
	2	1294.4530927801359	517.7812371120544	3	0.0	253.36721330595734	0.0	0.0	5396.589167868446;	%5 biomass
	2	0.0	0.0	3	0.0	0.42529944573960743	0.0	0.0	7085.78580603861;	%5 onwind
	2	0.0	0.0	3	0.0	0.3441067026914279	0.0	0.0	1934.5417707820682;	%5 solar
	2	133306.80450020512	53322.72180008205	3	0.0	680.1367576541077	0.0	0.0	43047.67045072198;	%6 CCGT
	2	188299.74316028284	75319.89726411313	3	0.0	960.7129753075653	0.0	0.0	19603.132567304394;	%6 OCGT
	2	1294.434965288245	517.7739861152979	3	0.0	253.36366515722153	0.0	0.0	5396.589167868446;	%6 biomass
	2	0.0	0.0	3	0.0	0.43126346650498476	0.0	0.0	7085.78580603861;	%6 onwind
	2	0.0	0.0	3	0.0	0.34342898895164936	0.0	0.0	1934.541770782068;	%6 solar
	2	135476.26918001496	54190.50767200598	3	0.0	691.2054550000762	0.0	0.0	43047.67045072198;	%7 CCGT
	2	189556.1329647723	75822.4531859089	3	0.0	967.1231273712871	0.0	0.0	19603.132567304394;	%7 OCGT
	2	1294.4519309635864	517.7807723854345	3	0.0	253.3669859000952	0.0	0.0	5396.589167868446;	%7 biomass
	2	0.0	0.0	3	0.0	0.42792709267635876	0.0	0.0	7085.78580603861;	%7 onwind
	2	0.0	0.0	3	0.0	0.33738894252321333	0.0	0.0	1934.5417707820682;	%7 solar
	2	145703.08067896435	58281.232271585744	3	0.0	743.3830646885937	0.0	0.0	43047.67045072198;	%8 CCGT
	2	190598.7295752882	76239.49183011528	3	0.0	972.4424978331028	0.0	0.0	19603.132567304394;	%8 OCGT
	2	1294.4620205251329	517.7848082100531	3	0.0	253.36896076044877	0.0	0.0	5396.589167868446;	%8 biomass
	2	0.0	0.0	3	0.0	0.43080802037477073	0.0	0.0	7085.78580603861;	%8 onwind
	2	0.0	0.0	3	0.0	0.339142701326348	0.0	0.0	1934.5417707820682;	%8 solar
	2	146733.38130083855	58693.352520335415	3	0.0	748.6397005144822	0.0	0.0	43047.67045072198;	%9 CCGT
	2	191316.0692204405	76526.4276881762	3	0.0	976.1023939818392	0.0	0.0	19603.132567304394;	%9 OCGT
	2	1294.411076654414	517.7644306617657	3	0.0	253.35898936277434	0.0	0.0	5396.589167868446;	%9 biomass
	2	0.0	0.0	3	0.0	0.4211632206701166	0.0	0.0	6712.888093376643;	%9 offwind
	2	0.0	0.0	3	0.0	0.42463168934036144	0.0	0.0	7085.78580603861;	%9 onwind
	2	0.0	0.0	3	0.0	0.3398536155007997	0.0	0.0	1934.5417707820682;	%9 solar
	2	0.0	0.0	3	0.0	0.1643574777513038	0.0	0.0	344.16309161315854;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.12326810831347784	0.0	0.0	344.16309161315854;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17906433862030915	0.0	0.0	344.16309161315854;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.13429825396523187	0.0	0.0	344.16309161315854;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18684403530930685	0.0	0.0	344.16309161315854;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.14013302648198014	0.0	0.0	344.16309161315854;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15603354962402194	0.0	0.0	344.16309161315854;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.11702516221801645	0.0	0.0	344.16309161315854;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18563467388490987	0.0	0.0	344.16309161315854;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.1392260054136824	0.0	0.0	344.16309161315854;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1785806819268018	0.0	0.0	344.16309161315854;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.13393551144510135	0.0	0.0	344.16309161315854;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16998510430562086	0.0	0.0	344.16309161315854;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.12748882822921564	0.0	0.0	344.16309161315854;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17254761969655857	0.0	0.0	344.16309161315854;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.17060553860900285	0.0	0.0	344.16309161315854;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.12795415395675214	0.0	0.0	344.16309161315854;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1661314004971328	0.0	0.0	344.16309161315854;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.12459855037284959	0.0	0.0	344.16309161315854;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1630778519377525	0.0	0.0	344.16309161315854;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.12230838895331439	0.0	0.0	344.16309161315854;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1561709299375128	0.0	0.0	344.16309161315854;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.1171281974531346	0.0	0.0	344.16309161315854;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1722516854311618	0.0	0.0	344.16309161315854;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.12918876407337135	0.0	0.0	344.16309161315854;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1581186307029683	0.0	0.0	344.16309161315854;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.11858897302722624	0.0	0.0	344.16309161315854;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15906990482249397	0.0	0.0	344.16309161315854;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.11930242861687049	0.0	0.0	344.16309161315854;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1843417965032332	0.0	0.0	344.16309161315854;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.1382563473774249	0.0	0.0	344.16309161315854;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1687506071389532	0.0	0.0	344.16309161315854;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.1265629553542149	0.0	0.0	344.16309161315854;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17538377389352428	0.0	0.0	344.16309161315854;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.13153783042014322	0.0	0.0	344.16309161315854;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1778037186734356	0.0	0.0	344.16309161315854;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.13335278900507672	0.0	0.0	344.16309161315854;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1578686862544206	0.0	0.0	344.16309161315854;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.18375747505802867	0.0	0.0	344.16309161315854;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.1378181062935215	0.0	0.0	344.16309161315854;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_003_009
%MPC_003_009	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 3 	Weight: 16
%	Time step: 9

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 16;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	3145.414835691814	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	11401.233022398972	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6820.399530994328	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	6572.939107421038	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	8741.319210705456	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	12446.808645138419	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	4862.450748268951	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	6047.798738379346	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6869.843597222378	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	7201.869764100594	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.109911892857827	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	27.42094051286295	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.855614413482354	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.733760871103	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.9984495715315442	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	23.163095691383436	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.716422664536425	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.713841235997505	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.242500240814312	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5448109206280978	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.046436942818275	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5122874051201736	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	21.379946640807628	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.933730860363334	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	24.517094780826763	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.668652316312959	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.61094324432382	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1658599570665196	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.7306035228985	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6724925285181764	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.01913922264071	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	24.55924703706284	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.7447026868320576	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	133439.35645694297	53375.74258277719	3	0.0	680.8130431476682	0.0	0.0	40515.45454185598;	%0 CCGT
	2	180011.9885457941	72004.79541831763	3	0.0	918.4285129887452	0.0	0.0	18450.00712216884;	%0 OCGT
	2	1218.3026239229262	487.3210495691705	3	0.0	238.46205204989747	0.0	0.0	5079.142746229125;	%0 biomass
	2	0.0	0.0	3	0.0	0.39533239048641416	0.0	0.0	5974.518063271553;	%0 offwind
	2	0.0	0.0	3	0.0	0.3932358209679051	0.0	0.0	6668.974876271633;	%0 onwind
	2	0.0	0.0	3	0.0	0.3146809894463291	0.0	0.0	1820.745196030182;	%0 solar
	2	134038.78265644438	53615.51306257775	3	0.0	683.8713400838999	0.0	0.0	40515.45454185598;	%1 CCGT
	2	172400.8261724471	68960.33046897883	3	0.0	879.5960519002401	0.0	0.0	18450.00712216884;	%1 OCGT
	2	1218.276968888761	487.3107875555044	3	0.0	238.457030512578	0.0	0.0	5079.142746229125;	%1 biomass
	2	0.0	0.0	3	0.0	0.40023321374282905	0.0	0.0	6668.974876271633;	%1 onwind
	2	0.0	0.0	3	0.0	0.32209040548784107	0.0	0.0	1820.745196030182;	%1 solar
	2	127779.53681324693	51111.81472529877	3	0.0	651.9364123124842	0.0	0.0	40515.45454185598;	%2 CCGT
	2	177139.84314028162	70855.93725611265	3	0.0	903.7747098993959	0.0	0.0	18450.00712216884;	%2 OCGT
	2	1218.2669230970375	487.30676923881504	3	0.0	238.45506421942406	0.0	0.0	5079.142746229125;	%2 biomass
	2	0.0	0.0	3	0.0	0.4018082454184173	0.0	0.0	6668.974876271633;	%2 onwind
	2	0.0	0.0	3	0.0	0.31937278986607287	0.0	0.0	1820.745196030182;	%2 solar
	2	118318.22470960322	47327.28988384129	3	0.0	603.6644117836898	0.0	0.0	40515.45454185598;	%3 CCGT
	2	169968.88720400236	67987.55488160094	3	0.0	867.18820002042	0.0	0.0	18450.00712216884;	%3 OCGT
	2	1218.3016074872776	487.320642994911	3	0.0	238.46185309987814	0.0	0.0	5079.142746229125;	%3 biomass
	2	0.0	0.0	3	0.0	0.4023179202111556	0.0	0.0	6117.393598518037;	%3 offwind
	2	0.0	0.0	3	0.0	0.39511844215176456	0.0	0.0	6668.974876271634;	%3 onwind
	2	0.0	0.0	3	0.0	0.31760507040530533	0.0	0.0	1820.745196030182;	%3 solar
	2	126693.58302335865	50677.43320934346	3	0.0	646.3958317518297	0.0	0.0	40515.45454185598;	%4 CCGT
	2	174856.5819309103	69942.63277236412	3	0.0	892.1254180148484	0.0	0.0	18450.00712216884;	%4 OCGT
	2	1218.3029278319186	487.3211711327674	3	0.0	238.46211153492237	0.0	0.0	5079.142746229125;	%4 biomass
	2	0.0	0.0	3	0.0	0.39617261173970336	0.0	0.0	6668.974876271633;	%4 onwind
	2	0.0	0.0	3	0.0	0.31288641662664113	0.0	0.0	1820.745196030182;	%4 solar
	2	131910.9020201775	52764.360808071004	3	0.0	673.0148062253954	0.0	0.0	40515.45454185598;	%5 CCGT
	2	175667.2941625954	70266.91766503817	3	0.0	896.2617049112012	0.0	0.0	18450.00712216884;	%5 OCGT
	2	1218.3087932048338	487.32351728193356	3	0.0	238.46325958207748	0.0	0.0	5079.142746229125;	%5 biomass
	2	0.0	0.0	3	0.0	0.4002818312843364	0.0	0.0	6668.974876271633;	%5 onwind
	2	0.0	0.0	3	0.0	0.3238651319448733	0.0	0.0	1820.745196030182;	%5 solar
	2	125465.22776489894	50186.091105959575	3	0.0	640.128713086219	0.0	0.0	40515.45454185598;	%6 CCGT
	2	177223.2876802662	70889.31507210647	3	0.0	904.2004473482967	0.0	0.0	18450.00712216884;	%6 OCGT
	2	1218.2917320359952	487.31669281439804	3	0.0	238.4599201479732	0.0	0.0	5079.142746229125;	%6 biomass
	2	0.0	0.0	3	0.0	0.40589502729880916	0.0	0.0	6668.974876271633;	%6 onwind
	2	0.0	0.0	3	0.0	0.3232272837191994	0.0	0.0	1820.7451960301817;	%6 solar
	2	127507.0768753082	51002.83075012328	3	0.0	650.546310588307	0.0	0.0	40515.45454185598;	%7 CCGT
	2	178405.7722021386	71362.30888085545	3	0.0	910.2335316435643	0.0	0.0	18450.00712216884;	%7 OCGT
	2	1218.3076997304343	487.32307989217367	3	0.0	238.46304555303075	0.0	0.0	5079.142746229125;	%7 biomass
	2	0.0	0.0	3	0.0	0.40275491075422	0.0	0.0	6668.974876271633;	%7 onwind
	2	0.0	0.0	3	0.0	0.3175425341394949	0.0	0.0	1820.745196030182;	%7 solar
	2	137132.31122726056	54852.92449090423	3	0.0	699.6546491186764	0.0	0.0	40515.45454185598;	%8 CCGT
	2	179387.03960027124	71754.8158401085	3	0.0	915.2399979605674	0.0	0.0	18450.00712216884;	%8 OCGT
	2	1218.3171957883603	487.3268783153441	3	0.0	238.46490424512825	0.0	0.0	5079.142746229125;	%8 biomass
	2	0.0	0.0	3	0.0	0.40546637211743125	0.0	0.0	6668.974876271633;	%8 onwind
	2	0.0	0.0	3	0.0	0.31919313066009225	0.0	0.0	1820.745196030182;	%8 solar
	2	138102.00593020098	55240.80237208039	3	0.0	704.6020710724539	0.0	0.0	40515.45454185598;	%9 CCGT
	2	180062.1827957087	72024.87311828349	3	0.0	918.6846061005546	0.0	0.0	18450.00712216884;	%9 OCGT
	2	1218.2692486159192	487.3076994463677	3	0.0	238.4555194002582	0.0	0.0	5079.142746229125;	%9 biomass
	2	0.0	0.0	3	0.0	0.39638891357187445	0.0	0.0	6318.012323178017;	%9 offwind
	2	0.0	0.0	3	0.0	0.39965335467328134	0.0	0.0	6668.974876271633;	%9 onwind
	2	0.0	0.0	3	0.0	0.3198622263536938	0.0	0.0	1820.745196030182;	%9 solar
	2	0.0	0.0	3	0.0	0.15468939082475652	0.0	0.0	323.91820387120805;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.11601704311856739	0.0	0.0	323.91820387120805;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1685311422308792	0.0	0.0	323.91820387120805;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.1263983566731594	0.0	0.0	323.91820387120805;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17585320970287704	0.0	0.0	323.91820387120805;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.13188990727715777	0.0	0.0	323.91820387120805;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.14685510552849124	0.0	0.0	323.91820387120805;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.11014132914636843	0.0	0.0	323.91820387120805;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17471498718579753	0.0	0.0	323.91820387120805;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.13103624038934814	0.0	0.0	323.91820387120805;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16807593593110756	0.0	0.0	323.91820387120805;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.12605695194833066	0.0	0.0	323.91820387120805;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15998598052293728	0.0	0.0	323.91820387120805;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.11998948539220296	0.0	0.0	323.91820387120805;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16239775971440806	0.0	0.0	323.91820387120805;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.1605699186908262	0.0	0.0	323.91820387120805;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.12042743901811966	0.0	0.0	323.91820387120805;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15635896517377204	0.0	0.0	323.91820387120805;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.11726922388032904	0.0	0.0	323.91820387120805;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15348503711788472	0.0	0.0	323.91820387120805;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.11511377783841353	0.0	0.0	323.91820387120805;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.14698440464707088	0.0	0.0	323.91820387120805;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.11023830348530317	0.0	0.0	323.91820387120805;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1621192333469758	0.0	0.0	323.91820387120805;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.12158942501023186	0.0	0.0	323.91820387120805;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1488175347792643	0.0	0.0	323.91820387120805;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.11161315108444822	0.0	0.0	323.91820387120805;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1497128515976414	0.0	0.0	323.91820387120805;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.11228463869823105	0.0	0.0	323.91820387120805;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17349816141480773	0.0	0.0	323.91820387120805;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.1301236210611058	0.0	0.0	323.91820387120805;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15882410083666185	0.0	0.0	323.91820387120805;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.11911807562749638	0.0	0.0	323.91820387120805;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16506708131155226	0.0	0.0	323.91820387120805;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.1238003109836642	0.0	0.0	323.91820387120805;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16734467639852763	0.0	0.0	323.91820387120805;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.12550850729889573	0.0	0.0	323.91820387120805;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.14858229294533704	0.0	0.0	323.91820387120805;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.1729482118193211	0.0	0.0	323.91820387120805;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.12971115886449083	0.0	0.0	323.91820387120805;	%DE1 76 PHS (pump mode)
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

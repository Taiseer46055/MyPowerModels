function mpc = mpc_000_013
%MPC_000_013	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 0 	Weight: 32
%	Time step: 13

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 32;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1745.707487565122	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	9794.31183097527	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	7261.087520089951	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3904.774912063142	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	5037.251940058762	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	4540.000088417577	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	10868.517583735742	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	6143.626056716054	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6730.637075127314	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	9137.722993019468	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.73205320287026	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.335869453119345	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.435431178083495	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.725582722743841	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.945645836302734	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.207064144289328	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.1356062140371375	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.552302555525024	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.807297700728212	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.768941597112836	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.7567071177994835	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.528522838859401	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.618546454146522	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.3802460996206225	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.631211531610179	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.583156536135236	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.10977044529962	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.637783127795208	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.6838236868100855	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.582657144160658	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.136317898421945	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.988649246050449	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.488256522745069	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	266878.71291388595	106751.48516555438	3	0.0	1361.6260862953363	0.0	0.0	81030.90908371196;	%0 CCGT
	2	360023.9770915882	144009.59083663527	3	0.0	1836.8570259774904	0.0	0.0	36900.01424433768;	%0 OCGT
	2	2436.6052478458523	974.642099138341	3	0.0	476.92410409979493	0.0	0.0	10158.28549245825;	%0 biomass
	2	0.0	0.0	3	0.0	0.7906647809728283	0.0	0.0	11949.036126543106;	%0 offwind
	2	0.0	0.0	3	0.0	0.7864716419358102	0.0	0.0	13337.949752543265;	%0 onwind
	2	0.0	0.0	3	0.0	0.6293619788926582	0.0	0.0	3641.490392060364;	%0 solar
	2	268077.56531288876	107231.0261251555	3	0.0	1367.7426801677998	0.0	0.0	81030.90908371196;	%1 CCGT
	2	344801.6523448942	137920.66093795767	3	0.0	1759.1921038004803	0.0	0.0	36900.01424433768;	%1 OCGT
	2	2436.553937777522	974.6215751110088	3	0.0	476.914061025156	0.0	0.0	10158.28549245825;	%1 biomass
	2	0.0	0.0	3	0.0	0.8004664274856581	0.0	0.0	13337.949752543265;	%1 onwind
	2	0.0	0.0	3	0.0	0.6441808109756821	0.0	0.0	3641.490392060364;	%1 solar
	2	255559.07362649386	102223.62945059754	3	0.0	1303.8728246249684	0.0	0.0	81030.90908371196;	%2 CCGT
	2	354279.68628056324	141711.8745122253	3	0.0	1807.5494197987919	0.0	0.0	36900.01424433768;	%2 OCGT
	2	2436.533846194075	974.6135384776301	3	0.0	476.9101284388481	0.0	0.0	10158.28549245825;	%2 biomass
	2	0.0	0.0	3	0.0	0.8036164908368346	0.0	0.0	13337.949752543265;	%2 onwind
	2	0.0	0.0	3	0.0	0.6387455797321457	0.0	0.0	3641.490392060364;	%2 solar
	2	236636.44941920644	94654.57976768258	3	0.0	1207.3288235673797	0.0	0.0	81030.90908371196;	%3 CCGT
	2	339937.7744080047	135975.10976320188	3	0.0	1734.37640004084	0.0	0.0	36900.01424433768;	%3 OCGT
	2	2436.603214974555	974.641285989822	3	0.0	476.9237061997563	0.0	0.0	10158.28549245825;	%3 biomass
	2	0.0	0.0	3	0.0	0.8046358404223112	0.0	0.0	12234.787197036074;	%3 offwind
	2	0.0	0.0	3	0.0	0.7902368843035291	0.0	0.0	13337.949752543267;	%3 onwind
	2	0.0	0.0	3	0.0	0.6352101408106107	0.0	0.0	3641.490392060364;	%3 solar
	2	253387.1660467173	101354.86641868691	3	0.0	1292.7916635036595	0.0	0.0	81030.90908371196;	%4 CCGT
	2	349713.1638618206	139885.26554472823	3	0.0	1784.2508360296968	0.0	0.0	36900.01424433768;	%4 OCGT
	2	2436.6058556638372	974.6423422655348	3	0.0	476.92422306984474	0.0	0.0	10158.28549245825;	%4 biomass
	2	0.0	0.0	3	0.0	0.7923452234794067	0.0	0.0	13337.949752543265;	%4 onwind
	2	0.0	0.0	3	0.0	0.6257728332532823	0.0	0.0	3641.490392060364;	%4 solar
	2	263821.804040355	105528.72161614201	3	0.0	1346.0296124507909	0.0	0.0	81030.90908371196;	%5 CCGT
	2	351334.5883251908	140533.83533007634	3	0.0	1792.5234098224023	0.0	0.0	36900.01424433768;	%5 OCGT
	2	2436.6175864096676	974.6470345638671	3	0.0	476.92651916415497	0.0	0.0	10158.28549245825;	%5 biomass
	2	0.0	0.0	3	0.0	0.8005636625686728	0.0	0.0	13337.949752543265;	%5 onwind
	2	0.0	0.0	3	0.0	0.6477302638897466	0.0	0.0	3641.490392060364;	%5 solar
	2	250930.45552979788	100372.18221191915	3	0.0	1280.257426172438	0.0	0.0	81030.90908371196;	%6 CCGT
	2	354446.5753605324	141778.63014421295	3	0.0	1808.4008946965935	0.0	0.0	36900.01424433768;	%6 OCGT
	2	2436.5834640719904	974.6333856287961	3	0.0	476.9198402959464	0.0	0.0	10158.28549245825;	%6 biomass
	2	0.0	0.0	3	0.0	0.8117900545976183	0.0	0.0	13337.949752543265;	%6 onwind
	2	0.0	0.0	3	0.0	0.6464545674383988	0.0	0.0	3641.4903920603633;	%6 solar
	2	255014.1537506164	102005.66150024656	3	0.0	1301.092621176614	0.0	0.0	81030.90908371196;	%7 CCGT
	2	356811.5444042772	142724.6177617109	3	0.0	1820.4670632871287	0.0	0.0	36900.01424433768;	%7 OCGT
	2	2436.6153994608685	974.6461597843473	3	0.0	476.9260911060615	0.0	0.0	10158.28549245825;	%7 biomass
	2	0.0	0.0	3	0.0	0.80550982150844	0.0	0.0	13337.949752543265;	%7 onwind
	2	0.0	0.0	3	0.0	0.6350850682789898	0.0	0.0	3641.490392060364;	%7 solar
	2	274264.62245452113	109705.84898180846	3	0.0	1399.3092982373528	0.0	0.0	81030.90908371196;	%8 CCGT
	2	358774.0792005425	143509.631680217	3	0.0	1830.4799959211348	0.0	0.0	36900.01424433768;	%8 OCGT
	2	2436.6343915767206	974.6537566306882	3	0.0	476.9298084902565	0.0	0.0	10158.28549245825;	%8 biomass
	2	0.0	0.0	3	0.0	0.8109327442348625	0.0	0.0	13337.949752543265;	%8 onwind
	2	0.0	0.0	3	0.0	0.6383862613201845	0.0	0.0	3641.490392060364;	%8 solar
	2	276204.01186040195	110481.60474416078	3	0.0	1409.2041421449078	0.0	0.0	81030.90908371196;	%9 CCGT
	2	360124.3655914174	144049.74623656698	3	0.0	1837.3692122011091	0.0	0.0	36900.01424433768;	%9 OCGT
	2	2436.5384972318384	974.6153988927354	3	0.0	476.9110388005164	0.0	0.0	10158.28549245825;	%9 biomass
	2	0.0	0.0	3	0.0	0.7927778271437489	0.0	0.0	12636.024646356034;	%9 offwind
	2	0.0	0.0	3	0.0	0.7993067093465627	0.0	0.0	13337.949752543265;	%9 onwind
	2	0.0	0.0	3	0.0	0.6397244527073876	0.0	0.0	3641.490392060364;	%9 solar
	2	0.0	0.0	3	0.0	0.30937878164951305	0.0	0.0	647.8364077424161;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.23203408623713478	0.0	0.0	647.8364077424161;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3370622844617584	0.0	0.0	647.8364077424161;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.2527967133463188	0.0	0.0	647.8364077424161;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3517064194057541	0.0	0.0	647.8364077424161;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.26377981455431554	0.0	0.0	647.8364077424161;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2937102110569825	0.0	0.0	647.8364077424161;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.22028265829273685	0.0	0.0	647.8364077424161;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34942997437159506	0.0	0.0	647.8364077424161;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.26207248077869627	0.0	0.0	647.8364077424161;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3361518718622151	0.0	0.0	647.8364077424161;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.25211390389666133	0.0	0.0	647.8364077424161;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31997196104587455	0.0	0.0	647.8364077424161;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.23997897078440592	0.0	0.0	647.8364077424161;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3247955194288161	0.0	0.0	647.8364077424161;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.3211398373816524	0.0	0.0	647.8364077424161;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.2408548780362393	0.0	0.0	647.8364077424161;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3127179303475441	0.0	0.0	647.8364077424161;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.23453844776065808	0.0	0.0	647.8364077424161;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30697007423576944	0.0	0.0	647.8364077424161;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.23022755567682707	0.0	0.0	647.8364077424161;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29396880929414176	0.0	0.0	647.8364077424161;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.22047660697060634	0.0	0.0	647.8364077424161;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3242384666939516	0.0	0.0	647.8364077424161;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.24317885002046372	0.0	0.0	647.8364077424161;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2976350695585286	0.0	0.0	647.8364077424161;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.22322630216889644	0.0	0.0	647.8364077424161;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2994257031952828	0.0	0.0	647.8364077424161;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.2245692773964621	0.0	0.0	647.8364077424161;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34699632282961546	0.0	0.0	647.8364077424161;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.2602472421222116	0.0	0.0	647.8364077424161;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3176482016733237	0.0	0.0	647.8364077424161;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.23823615125499276	0.0	0.0	647.8364077424161;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3301341626231045	0.0	0.0	647.8364077424161;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.2476006219673284	0.0	0.0	647.8364077424161;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33468935279705525	0.0	0.0	647.8364077424161;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.25101701459779147	0.0	0.0	647.8364077424161;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29716458589067407	0.0	0.0	647.8364077424161;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.3458964236386422	0.0	0.0	647.8364077424161;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.25942231772898167	0.0	0.0	647.8364077424161;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_011_008
%MPC_011_008	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 11 	Weight: 11
%	Time step: 8

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 11;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	2965.420622716881	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	10016.575215751145	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	7351.371952148461	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	6076.804743236598	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	12186.458583641059	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	14784.41678693841	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3604.427904861039	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	6063.715009568495	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	8289.476797121137	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	6570.819415413156	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.676875448233048	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	36.07672759458159	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.3681616225096498	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	35.550361145873396	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.058577549846693	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	35.80466464621714	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.393884683951013	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.147849634501252	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	37.23178524557011	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.5086246693693426	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	30.362241850247585	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.766241299698505	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	33.19883432593393	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6487972957525576	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	37.73538597976095	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.331143579431871	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	36.2499520684895	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2000533733357837	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	35.59837557658188	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0659161679098155	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.191360544749529	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	37.926387945668246	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.727343675077886	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	91739.5575641483	36695.823025659316	3	0.0	468.05896716402185	0.0	0.0	27854.374997525985;	%0 CCGT
	2	123758.24212523345	49503.29685009337	3	0.0	631.4196026797623	0.0	0.0	12684.379896491078;	%0 OCGT
	2	837.5830539470118	335.0332215788047	3	0.0	163.9426607843045	0.0	0.0	3491.9106380325234;	%0 biomass
	2	0.0	0.0	3	0.0	0.27179101845940973	0.0	0.0	4107.481168499193;	%0 offwind
	2	0.0	0.0	3	0.0	0.27034962691543474	0.0	0.0	4584.920227436747;	%0 onwind
	2	0.0	0.0	3	0.0	0.21634318024435126	0.0	0.0	1251.76232227075;	%0 solar
	2	92151.66307630551	36860.66523052221	3	0.0	470.1615463076812	0.0	0.0	27854.374997525985;	%1 CCGT
	2	118525.56799355737	47410.22719742295	3	0.0	604.7222856814151	0.0	0.0	12684.379896491078;	%1 OCGT
	2	837.5654161110233	335.0261664444093	3	0.0	163.9392084773974	0.0	0.0	3491.9106380325234;	%1 biomass
	2	0.0	0.0	3	0.0	0.275160334448195	0.0	0.0	4584.920227436747;	%1 onwind
	2	0.0	0.0	3	0.0	0.22143715377289075	0.0	0.0	1251.76232227075;	%1 solar
	2	87848.43155910727	35139.37262364291	3	0.0	448.2062834648329	0.0	0.0	27854.374997525985;	%2 CCGT
	2	121783.6421589436	48713.45686357745	3	0.0	621.3451130558346	0.0	0.0	12684.379896491078;	%2 OCGT
	2	837.5585096292133	335.0234038516853	3	0.0	163.93785665085403	0.0	0.0	3491.9106380325234;	%2 biomass
	2	0.0	0.0	3	0.0	0.2762431687251619	0.0	0.0	4584.920227436747;	%2 onwind
	2	0.0	0.0	3	0.0	0.2195687930329251	0.0	0.0	1251.76232227075;	%2 solar
	2	81343.77948785221	32537.511795140887	3	0.0	415.01928310128676	0.0	0.0	27854.374997525985;	%3 CCGT
	2	116853.60995275162	46741.443981100645	3	0.0	596.1918875140387	0.0	0.0	12684.379896491078;	%3 OCGT
	2	837.5823551475033	335.0329420590013	3	0.0	163.94252400616622	0.0	0.0	3491.9106380325234;	%3 biomass
	2	0.0	0.0	3	0.0	0.27659357014516944	0.0	0.0	4205.70809898115;	%3 offwind
	2	0.0	0.0	3	0.0	0.27164392897933815	0.0	0.0	4584.920227436748;	%3 onwind
	2	0.0	0.0	3	0.0	0.21835348590364742	0.0	0.0	1251.76232227075;	%3 solar
	2	87101.83832855907	34840.735331423624	3	0.0	444.39713432938294	0.0	0.0	27854.374997525985;	%4 CCGT
	2	120213.90007750083	48085.56003100033	3	0.0	613.3362248852083	0.0	0.0	12684.379896491078;	%4 OCGT
	2	837.583262884444	335.0333051537776	3	0.0	163.94270168025912	0.0	0.0	3491.9106380325234;	%4 biomass
	2	0.0	0.0	3	0.0	0.2723686705710461	0.0	0.0	4584.920227436747;	%4 onwind
	2	0.0	0.0	3	0.0	0.21510941143081577	0.0	0.0	1251.76232227075;	%4 solar
	2	90688.74513887204	36275.498055548815	3	0.0	462.69767927995935	0.0	0.0	27854.374997525985;	%5 CCGT
	2	120771.26473678435	48308.505894713744	3	0.0	616.1799221264508	0.0	0.0	12684.379896491078;	%5 OCGT
	2	837.5872953283233	335.03491813132933	3	0.0	163.94349096267828	0.0	0.0	3491.9106380325234;	%5 biomass
	2	0.0	0.0	3	0.0	0.27519375900798126	0.0	0.0	4584.920227436747;	%5 onwind
	2	0.0	0.0	3	0.0	0.2226572782121004	0.0	0.0	1251.76232227075;	%5 solar
	2	86257.34408836802	34502.93763534721	3	0.0	440.08849024677556	0.0	0.0	27854.374997525985;	%6 CCGT
	2	121841.010280183	48736.4041120732	3	0.0	621.637807551954	0.0	0.0	12684.379896491078;	%6 OCGT
	2	837.5755657747467	335.03022630989864	3	0.0	163.94119510173158	0.0	0.0	3491.9106380325234;	%6 biomass
	2	0.0	0.0	3	0.0	0.2790528312679313	0.0	0.0	4584.920227436747;	%6 onwind
	2	0.0	0.0	3	0.0	0.2222187575569496	0.0	0.0	1251.7623222707498;	%6 solar
	2	87661.11535177438	35064.446140709755	3	0.0	447.2505885294611	0.0	0.0	27854.374997525985;	%7 CCGT
	2	122653.9683889703	49061.58735558812	3	0.0	625.7855530049504	0.0	0.0	12684.379896491078;	%7 OCGT
	2	837.5865435646735	335.0346174258694	3	0.0	163.94334381770864	0.0	0.0	3491.9106380325234;	%7 biomass
	2	0.0	0.0	3	0.0	0.27689400114352625	0.0	0.0	4584.920227436747;	%7 onwind
	2	0.0	0.0	3	0.0	0.21831049222090274	0.0	0.0	1251.76232227075;	%7 solar
	2	94278.46396874163	37711.38558749666	3	0.0	481.01257126909	0.0	0.0	27854.374997525985;	%8 CCGT
	2	123328.58972518648	49331.43589007459	3	0.0	629.2274985978901	0.0	0.0	12684.379896491078;	%8 OCGT
	2	837.5930721044978	335.0372288417991	3	0.0	163.94462166852568	0.0	0.0	3491.9106380325234;	%8 biomass
	2	0.0	0.0	3	0.0	0.278758130830734	0.0	0.0	4584.920227436747;	%8 onwind
	2	0.0	0.0	3	0.0	0.21944527732881342	0.0	0.0	1251.76232227075;	%8 solar
	2	94945.12907701317	37978.05163080527	3	0.0	484.41392386231206	0.0	0.0	27854.374997525985;	%9 CCGT
	2	123792.75067204973	49517.1002688199	3	0.0	631.5956666941313	0.0	0.0	12684.379896491078;	%9 OCGT
	2	837.5601084234445	335.02404336937775	3	0.0	163.93816958767752	0.0	0.0	3491.9106380325234;	%9 biomass
	2	0.0	0.0	3	0.0	0.2725173780806637	0.0	0.0	4343.633472184887;	%9 offwind
	2	0.0	0.0	3	0.0	0.2747616813378809	0.0	0.0	4584.920227436747;	%9 onwind
	2	0.0	0.0	3	0.0	0.2199052806181645	0.0	0.0	1251.76232227075;	%9 solar
	2	0.0	0.0	3	0.0	0.10634895619202012	0.0	0.0	222.69376516145553;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.07976171714401509	0.0	0.0	222.69376516145553;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11586516028372945	0.0	0.0	222.69376516145553;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.08689887021279709	0.0	0.0	222.69376516145553;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12089908167072796	0.0	0.0	222.69376516145553;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.09067431125304598	0.0	0.0	222.69376516145553;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10096288505083773	0.0	0.0	222.69376516145553;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.0757221637881283	0.0	0.0	222.69376516145553;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1201165536902358	0.0	0.0	222.69376516145553;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.09008741526767686	0.0	0.0	222.69376516145553;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11555220595263645	0.0	0.0	222.69376516145553;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.08666415446447734	0.0	0.0	222.69376516145553;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10999036160951937	0.0	0.0	222.69376516145553;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.08249277120713953	0.0	0.0	222.69376516145553;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11164845980365554	0.0	0.0	222.69376516145553;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.11039181909994301	0.0	0.0	222.69376516145553;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.08279386432495726	0.0	0.0	222.69376516145553;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10749678855696827	0.0	0.0	222.69376516145553;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.08062259141772621	0.0	0.0	222.69376516145553;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10552096301854574	0.0	0.0	222.69376516145553;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.07914072226390931	0.0	0.0	222.69376516145553;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10105177819486123	0.0	0.0	222.69376516145553;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.07578883364614593	0.0	0.0	222.69376516145553;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11145697292604587	0.0	0.0	222.69376516145553;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.0835927296945344	0.0	0.0	222.69376516145553;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10231205516074421	0.0	0.0	222.69376516145553;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.07673404137055816	0.0	0.0	222.69376516145553;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10292758547337845	0.0	0.0	222.69376516145553;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.07719568910503384	0.0	0.0	222.69376516145553;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11927998597268032	0.0	0.0	222.69376516145553;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.08945998947951024	0.0	0.0	222.69376516145553;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10919156932520502	0.0	0.0	222.69376516145553;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.08189367699390376	0.0	0.0	222.69376516145553;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11348361840169217	0.0	0.0	222.69376516145553;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.08511271380126913	0.0	0.0	222.69376516145553;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11504946502398775	0.0	0.0	222.69376516145553;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.08628709876799082	0.0	0.0	222.69376516145553;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.10215032639991921	0.0	0.0	222.69376516145553;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.11890189562578325	0.0	0.0	222.69376516145553;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.08917642171933744	0.0	0.0	222.69376516145553;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_014_013
%MPC_014_013	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 14 	Weight: 22
%	Time step: 13

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 22;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1583.6206576621305	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5838.4741854064	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1161.9893190009036	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-121.87376539275272	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	3331.707377961266	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	7570.226757228691	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	302.1871486034454	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3093.6664357966197	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2699.883467017427	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	4208.482532154018	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	-447.8477113391523	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	4440.443239912751	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6060.862919032055	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	5019.69606962486	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	4844.512261636997	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	4608.426769220448	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	6316.596473604902	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5028.794297855558	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5370.700838736096	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4309.01586057852	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.524400116473826	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	23.484262348080225	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.36521339828071	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.254613739392623	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.2955940549787672	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.948003778923981	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	21.72530525418379	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.173957977351286	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.265073238893228	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	34.381597248567424	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.430339684931716	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.996187583362037	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.39309474475815914	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.3631244713539745	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.793732234521515	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	35.23200247980568	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.782873314477948	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.264706644958313	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	28.22606036721543	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1068721099658316	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.3163734739939064	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.321178210145712	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5212408209116873	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.914019722254542	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.0210807947796	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.333462421096128	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.306081506185024	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3587882699075977	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.137447887227023	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.0353303280611024	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.970330217593327	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	22.755853331738283	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2471099590425263	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.08141531711466476	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.035396198304849	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.503598680942409	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.60195882613325	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.904753430582662	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.897880232537277	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0582745684378887	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.63232656771274	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2207023607650433	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.822378477312139	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.3309672293759998	0.0	0.0	132	132	890	28.0	3;	%9 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1776.8	266.52	3.5	1	1	1	1776.8	7;	%DE1 0 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-266.52	-1776.8	3.5	1	1	1	-1776.8	8;	%DE1 0 PHS (pump mode)
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	1294.0	194.1	3.5	1	1	1	1294.0	7;	%DE1 12 PHS
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	-194.1	-1294.0	3.5	1	1	1	-1294.0	8;	%DE1 12 PHS (pump mode)
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
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	167.1	25.064999999999998	3.5	1	1	1	167.1	7;	%DE1 27 PHS
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	-25.064999999999998	-167.1	3.5	1	1	1	-167.1	8;	%DE1 27 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	109.0	16.349999999999998	3.5	1	1	1	109.0	7;	%DE1 31 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-16.349999999999998	-109.0	3.5	1	1	1	-109.0	8;	%DE1 31 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	480.0	72.0	3.5	1	1	1	480.0	7;	%DE1 33 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-72.0	-480.0	3.5	1	1	1	-480.0	8;	%DE1 33 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	92.0	13.799999999999999	3.5	1	1	1	92.0	7;	%DE1 36 PHS
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.799999999999999	-92.0	3.5	1	1	1	-92.0	8;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	278.0	41.699999999999996	3.5	1	1	1	278.0	7;	%DE1 43 PHS
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	-41.699999999999996	-278.0	3.5	1	1	1	-278.0	8;	%DE1 43 PHS (pump mode)
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	164.0	24.599999999999998	3.5	1	1	1	164.0	7;	%DE1 49 PHS
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.599999999999998	-164.0	3.5	1	1	1	-164.0	8;	%DE1 49 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.7	11.955	3.5	1	1	1	79.7	7;	%DE1 50 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-11.955	-79.7	3.5	1	1	1	-79.7	8;	%DE1 50 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	399.8	59.97	3.5	1	1	1	399.8	7;	%DE1 52 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-59.97	-399.8	3.5	1	1	1	-399.8	8;	%DE1 52 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	162.0	24.3	3.5	1	1	1	162.0	7;	%DE1 54 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.3	-162.0	3.5	1	1	1	-162.0	8;	%DE1 54 PHS (pump mode)
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	90.0	13.5	3.5	1	1	1	90.0	7;	%DE1 61 PHS
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.5	-90.0	3.5	1	1	1	-90.0	8;	%DE1 61 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	1060.0	159.0	3.5	1	1	1	1060.0	7;	%DE1 7 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-159.0	-1060.0	3.5	1	1	1	-1060.0	8;	%DE1 7 PHS (pump mode)
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	45.5	6.825	3.5	1	1	1	45.5	6;	%DE1 71 hydro
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	18.0	3.5	1	1	1	120.0	7;	%DE1 76 PHS
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	-18.0	-120.0	3.5	1	1	1	-120.0	8;	%DE1 76 PHS (pump mode)
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0	investment
mpc.gencost = [
	2	1675.1973508742078	670.0789403496831	3	0.0	327.8914368514793	0.0	0.0	6983.821276065048;	%0 biomass
	2	0.0	0.0	3	0.0	0.5428038513305424	0.0	0.0	8088.366798979818;	%0 offwind
	2	0.0	0.0	3	0.0	0.5301446090267016	0.0	0.0	9169.840454873494;	%0 onwind
	2	0.0	0.0	3	0.0	0.41929493554560837	0.0	0.0	2503.5246445415;	%0 solar
	2	195423.86142924527	78169.5445716981	3	0.0	997.0605174961491	0.0	0.0	55708.74999505197;	%1 CCGT
	2	238545.94876672293	95418.37950668916	3	0.0	1217.0711671771576	0.0	0.0	25368.759792982157;	%1 OCGT
	2	1675.1304300829254	670.0521720331701	3	0.0	327.87833824289004	0.0	0.0	6983.821276065047;	%1 biomass
	2	0.0	0.0	3	0.0	0.5407603078217602	0.0	0.0	9169.840454873494;	%1 onwind
	2	0.0	0.0	3	0.0	0.4386749100713146	0.0	0.0	2503.5246445415;	%1 solar
	2	171232.09741394283	68492.83896557713	3	0.0	873.6331500711367	0.0	0.0	55708.74999505197;	%10 CCGT
	2	248384.5636481726	99353.82545926904	3	0.0	1267.2681818784317	0.0	0.0	25368.759792982157;	%10 OCGT
	2	1675.161245736134	670.0644982944536	3	0.0	327.88436988376077	0.0	0.0	6983.821276065047;	%10 biomass
	2	0.0	0.0	3	0.0	0.5432479479108692	0.0	0.0	7985.071707084331;	%10 offwind
	2	0.0	0.0	3	0.0	0.5455936039296327	0.0	0.0	9169.840454873496;	%10 onwind
	2	0.0	0.0	3	0.0	0.4360901162006844	0.0	0.0	2503.5246445415;	%10 solar
	2	159776.76276034096	63910.70510413638	3	0.0	815.1875651037803	0.0	0.0	55708.74999505197;	%11 CCGT
	2	237188.70738477254	94875.48295390903	3	0.0	1210.1464662488393	0.0	0.0	25368.759792982157;	%11 OCGT
	2	1675.1982934867692	670.0793173947077	3	0.0	327.8916213518828	0.0	0.0	6983.821276065047;	%11 biomass
	2	0.0	0.0	3	0.0	0.5569128597654622	0.0	0.0	7939.363060933662;	%11 offwind
	2	0.0	0.0	3	0.0	0.5303498265504576	0.0	0.0	9169.840454873494;	%11 onwind
	2	0.0	0.0	3	0.0	0.44209885630589796	0.0	0.0	2503.5246445415;	%11 solar
	2	235910.47989335994	94364.19195734397	3	0.0	1203.6248974151015	0.0	0.0	25368.759792982157;	%12 OCGT
	2	1675.1664552562092	670.0665821024837	3	0.0	327.8853895588587	0.0	0.0	6983.821276065047;	%12 biomass
	2	0.0	0.0	3	0.0	0.5503509640900184	0.0	0.0	9169.840454873494;	%12 onwind
	2	0.0	0.0	3	0.0	0.44022235933731024	0.0	0.0	2503.5246445415;	%12 solar
	2	177360.7192619493	70944.28770477972	3	0.0	904.9016288874964	0.0	0.0	55708.74999505197;	%13 CCGT
	2	246716.5567159531	98686.62268638123	3	0.0	1258.7579424283322	0.0	0.0	25368.759792982157;	%13 OCGT
	2	1675.1510782490877	670.0604312996352	3	0.0	327.8823797708138	0.0	0.0	6983.821276065047;	%13 biomass
	2	0.0	0.0	3	0.0	0.5594432523895452	0.0	0.0	9169.840454873494;	%13 onwind
	2	0.0	0.0	3	0.0	0.4453601277180325	0.0	0.0	2503.5246445415;	%13 solar
	2	182769.31204267306	73107.72481706922	3	0.0	932.4964900136381	0.0	0.0	55708.74999505197;	%14 CCGT
	2	246539.92288695704	98615.96915478281	3	0.0	1257.8567494232502	0.0	0.0	25368.759792982157;	%14 OCGT
	2	1675.0824162067079	670.0329664826832	3	0.0	327.8689403418884	0.0	0.0	6983.821276065047;	%14 biomass
	2	0.0	0.0	3	0.0	0.5439436278198639	0.0	0.0	8361.556536460052;	%14 offwind
	2	0.0	0.0	3	0.0	0.538560055047701	0.0	0.0	9169.840454873496;	%14 onwind
	2	0.0	0.0	3	0.0	0.4305370119127677	0.0	0.0	2503.5246445415;	%14 solar
	2	164142.95708338617	65657.18283335447	3	0.0	837.4640667519701	0.0	0.0	55708.74999505197;	%15 CCGT
	2	230225.73242623388	92090.29297049357	3	0.0	1174.621083807316	0.0	0.0	25368.759792982157;	%15 OCGT
	2	1675.1701747485058	670.0680698994023	3	0.0	327.8861175863194	0.0	0.0	6983.821276065047;	%15 biomass
	2	0.0	0.0	3	0.0	0.5507033273069233	0.0	0.0	8726.118289314727;	%15 offwind
	2	0.0	0.0	3	0.0	0.5451320664816096	0.0	0.0	9169.840454873496;	%15 onwind
	2	0.0	0.0	3	0.0	0.42546731928589265	0.0	0.0	2503.5246445415;	%15 solar
	2	171820.7166419517	68728.28665678068	3	0.0	876.6363093977126	0.0	0.0	55708.74999505197;	%16 CCGT
	2	237734.9993076775	95093.999723071	3	0.0	1212.93366993713	0.0	0.0	25368.759792982157;	%16 OCGT
	2	1675.1298086374322	670.0519234549729	3	0.0	327.87821660548684	0.0	0.0	6983.821276065048;	%16 biomass
	2	0.0	0.0	3	0.0	0.5442087902394104	0.0	0.0	9169.840454873494;	%16 onwind
	2	0.0	0.0	3	0.0	0.4376309853546966	0.0	0.0	2503.5246445415005;	%16 solar
	2	176654.0943213077	70661.63772852308	3	0.0	901.2963995985086	0.0	0.0	55708.74999505197;	%17 CCGT
	2	231753.47509642644	92701.39003857058	3	0.0	1182.4156892674816	0.0	0.0	25368.759792982157;	%17 OCGT
	2	1675.1487254576257	670.0594901830502	3	0.0	327.8819192518351	0.0	0.0	6983.821276065047;	%17 biomass
	2	0.0	0.0	3	0.0	0.5608265074371094	0.0	0.0	9169.840454873496;	%17 onwind
	2	0.0	0.0	3	0.0	0.44844781224308683	0.0	0.0	2503.5246445415005;	%17 solar
	2	185173.3043502758	74069.32174011033	3	0.0	944.7617568891623	0.0	0.0	55708.74999505197;	%18 CCGT
	2	233557.12130807556	93422.84852323022	3	0.0	1191.6179658575281	0.0	0.0	25368.759792982157;	%18 OCGT
	2	1675.1790454923396	670.0716181969358	3	0.0	327.88785388380103	0.0	0.0	6983.821276065047;	%18 biomass
	2	0.0	0.0	3	0.0	0.5480110962909743	0.0	0.0	9169.840454873494;	%18 onwind
	2	0.0	0.0	3	0.0	0.43998755306671467	0.0	0.0	2503.5246445415;	%18 solar
	2	173509.96929294264	69403.98771717705	3	0.0	885.2549453721563	0.0	0.0	55708.74999505197;	%19 CCGT
	2	243103.1488708119	97241.25954832477	3	0.0	1240.3221881163872	0.0	0.0	25368.759792982157;	%19 OCGT
	2	1675.1404461074374	670.0561784429749	3	0.0	327.88029870961776	0.0	0.0	6983.821276065047;	%19 biomass
	2	0.0	0.0	3	0.0	0.559108788952341	0.0	0.0	9169.840454873494;	%19 onwind
	2	0.0	0.0	3	0.0	0.43734160415765205	0.0	0.0	2503.5246445415;	%19 solar
	2	185758.80250070174	74303.52100028071	3	0.0	947.7489923505192	0.0	0.0	55708.74999505197;	%2 CCGT
	2	245503.42244074392	98201.36897629757	3	0.0	1252.56848184053	0.0	0.0	25368.759792982157;	%2 OCGT
	2	1675.132253686687	670.0529014746747	3	0.0	327.8786951823619	0.0	0.0	6983.821276065047;	%2 biomass
	2	0.0	0.0	3	0.0	0.5509980178638266	0.0	0.0	9169.840454873494;	%2 onwind
	2	0.0	0.0	3	0.0	0.4387464030011597	0.0	0.0	2503.5246445415;	%2 solar
	2	202574.27421871357	81029.70968748543	3	0.0	1033.5422154015998	0.0	0.0	55708.74999505197;	%3 CCGT
	2	251987.29172680006	100794.91669072003	3	0.0	1285.6494475857144	0.0	0.0	25368.759792982157;	%3 OCGT
	2	1675.0833201039993	670.0333280415997	3	0.0	327.8691172644352	0.0	0.0	6983.821276065047;	%3 biomass
	2	0.0	0.0	3	0.0	0.5450347561613273	0.0	0.0	8687.266944369774;	%3 offwind
	2	0.0	0.0	3	0.0	0.5607121753916999	0.0	0.0	9169.840454873496;	%3 onwind
	2	0.0	0.0	3	0.0	0.4393433526788577	0.0	0.0	2503.5246445415;	%3 solar
	2	176586.6366722846	70634.65466891385	3	0.0	900.9522279198194	0.0	0.0	55708.74999505197;	%4 CCGT
	2	241497.29674112977	96598.9186964519	3	0.0	1232.1290650057642	0.0	0.0	25368.759792982157;	%4 OCGT
	2	1675.182813615864	670.0731254463456	3	0.0	327.8885914299988	0.0	0.0	6983.821276065047;	%4 biomass
	2	0.0	0.0	3	0.0	0.5484557677877067	0.0	0.0	9169.840454873494;	%4 onwind
	2	0.0	0.0	3	0.0	0.42997985040275555	0.0	0.0	2503.5246445415;	%4 solar
	2	181377.49027774407	72550.99611109763	3	0.0	925.3953585599187	0.0	0.0	55708.74999505197;	%5 CCGT
	2	237145.77156480442	94858.30862592177	3	0.0	1209.9274059428794	0.0	0.0	25368.759792982157;	%5 OCGT
	2	1675.1605255056086	670.0642102022434	3	0.0	327.88422891086486	0.0	0.0	6983.821276065047;	%5 biomass
	2	0.0	0.0	3	0.0	0.5417880023153013	0.0	0.0	9169.840454873494;	%5 onwind
	2	0.0	0.0	3	0.0	0.4403741200742285	0.0	0.0	2503.5246445415;	%5 solar
	2	177838.84890597698	71135.53956239078	3	0.0	907.3410658468213	0.0	0.0	55708.74999505197;	%6 CCGT
	2	245615.96003223164	98246.38401289265	3	0.0	1253.1426532256717	0.0	0.0	25368.759792982157;	%6 OCGT
	2	1675.1883030443648	670.0753212177459	3	0.0	327.88966589241824	0.0	0.0	6983.821276065047;	%6 biomass
	2	0.0	0.0	3	0.0	0.5502132087275907	0.0	0.0	9169.840454873494;	%6 onwind
	2	0.0	0.0	3	0.0	0.44921485188695087	0.0	0.0	2503.5246445415;	%6 solar
	2	168565.48246727663	67426.19298691065	3	0.0	860.0279717718195	0.0	0.0	55708.74999505197;	%7 CCGT
	2	244565.90055183362	97826.36022073343	3	0.0	1247.7852068971101	0.0	0.0	25368.759792982157;	%7 OCGT
	2	1675.1865091444085	670.0746036577634	3	0.0	327.88931476696195	0.0	0.0	6983.821276065047;	%7 biomass
	2	0.0	0.0	3	0.0	0.5562314018244966	0.0	0.0	9169.840454873494;	%7 onwind
	2	0.0	0.0	3	0.0	0.43890228993005925	0.0	0.0	2503.5246445415005;	%7 solar
	2	190587.1020898078	76234.84083592311	3	0.0	972.3831739275906	0.0	0.0	55708.74999505197;	%8 CCGT
	2	249932.1939859473	99972.87759437892	3	0.0	1275.1642550303434	0.0	0.0	25368.759792982157;	%8 OCGT
	2	1675.189693567323	670.0758774269292	3	0.0	327.8899380636764	0.0	0.0	6983.821276065047;	%8 biomass
	2	0.0	0.0	3	0.0	0.5622688443467148	0.0	0.0	9169.840454873494;	%8 onwind
	2	0.0	0.0	3	0.0	0.4383420554530829	0.0	0.0	2503.5246445415;	%8 solar
	2	182756.87234215654	73102.74893686263	3	0.0	932.4330221538598	0.0	0.0	55708.74999505197;	%9 CCGT
	2	246589.34103001253	98635.73641200502	3	0.0	1258.1088828061863	0.0	0.0	25368.759792982157;	%9 OCGT
	2	1675.123950079237	670.0495800316949	3	0.0	327.8770698921975	0.0	0.0	6983.821276065047;	%9 biomass
	2	0.0	0.0	3	0.0	0.5468261079312333	0.0	0.0	9169.840454873494;	%9 onwind
	2	0.0	0.0	3	0.0	0.4434508286613258	0.0	0.0	2503.5246445415;	%9 solar
	2	0.0	0.0	3	0.0	0.21269791238404023	0.0	0.0	445.38753032291106;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.15952343428803017	0.0	0.0	445.38753032291106;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2317303205674589	0.0	0.0	445.38753032291106;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.17379774042559418	0.0	0.0	445.38753032291106;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24179816334145593	0.0	0.0	445.38753032291106;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.18134862250609196	0.0	0.0	445.38753032291106;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20192577010167545	0.0	0.0	445.38753032291106;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.1514443275762566	0.0	0.0	445.38753032291106;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2402331073804716	0.0	0.0	445.38753032291106;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.1801748305353537	0.0	0.0	445.38753032291106;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2311044119052729	0.0	0.0	445.38753032291106;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.17332830892895468	0.0	0.0	445.38753032291106;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21998072321903875	0.0	0.0	445.38753032291106;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.16498554241427907	0.0	0.0	445.38753032291106;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22329691960731107	0.0	0.0	445.38753032291106;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.22078363819988603	0.0	0.0	445.38753032291106;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.16558772864991453	0.0	0.0	445.38753032291106;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21499357711393655	0.0	0.0	445.38753032291106;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.16124518283545242	0.0	0.0	445.38753032291106;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2110419260370915	0.0	0.0	445.38753032291106;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.15828144452781862	0.0	0.0	445.38753032291106;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20210355638972247	0.0	0.0	445.38753032291106;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.15157766729229186	0.0	0.0	445.38753032291106;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22291394585209173	0.0	0.0	445.38753032291106;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.1671854593890688	0.0	0.0	445.38753032291106;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20462411032148842	0.0	0.0	445.38753032291106;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.1534680827411163	0.0	0.0	445.38753032291106;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2058551709467569	0.0	0.0	445.38753032291106;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.15439137821006768	0.0	0.0	445.38753032291106;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23855997194536063	0.0	0.0	445.38753032291106;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.17891997895902048	0.0	0.0	445.38753032291106;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21838313865041004	0.0	0.0	445.38753032291106;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.16378735398780753	0.0	0.0	445.38753032291106;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22696723680338435	0.0	0.0	445.38753032291106;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.17022542760253825	0.0	0.0	445.38753032291106;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2300989300479755	0.0	0.0	445.38753032291106;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.17257419753598163	0.0	0.0	445.38753032291106;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20430065279983842	0.0	0.0	445.38753032291106;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.2378037912515665	0.0	0.0	445.38753032291106;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.17835284343867489	0.0	0.0	445.38753032291106;	%DE1 76 PHS (pump mode)
];
%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	15	0.004446514507151612	0.03646141895864322	3.1100887306422827e-06	983.1120383760949	983.1120383760949	983.1120383760949	0	0	1	-60	60;
	1	9	0.0066236548409098605	0.05431396969546085	4.632876884410344e-06	983.1120383760949	983.1120383760949	983.1120383760949	0	0	1	-60	60;
	2	14	0.00036342889364109445	0.002980116927856975	4.3562433442329304e-05	12869.830320559791	12869.830320559791	12869.830320559791	0	0	1	-60	60;
	2	18	0.0005977660653844128	0.004901681736152186	1.995841906294707e-05	6792.410446962112	6792.410446962112	6792.410446962112	0	0	1	-60	60;
	2	20	0.0007594657170519633	0.0062276188798261	2.0299890501867426e-05	6077.419873597678	6077.419873597678	6077.419873597678	0	0	1	-60	60;
	2	7	0.0010122549625824677	0.008300490693176236	2.7056737965948302e-05	6077.419873597678	6077.419873597678	6077.419873597678	0	0	1	-60	60;
	2	8	0.002010615542990133	0.016487047452519093	1.6782780884218432e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	2	10	0.0005241817991184883	0.004298290752771604	5.3598534430142926e-05	11886.718282183694	11886.718282183694	11886.718282183694	0	0	1	-60	60;
	11	15	0.0011444155077600761	0.009384207163632623	2.381518341854371e-05	5362.429300233245	5362.429300233245	5362.429300233245	0	0	1	-60	60;
	11	16	0.0014237338026427768	0.01167461718167077	1.1884028515801675e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	11	3	0.003092914165728497	0.025361896158973676	2.581682058416841e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	11	8	0.0021722130471396196	0.017812146986544878	1.8131649151467633e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	11	9	0.0008540008001274558	0.007002806561045137	2.8513672548599497e-05	6792.410446962112	6792.410446962112	6792.410446962112	0	0	1	-60	60;
	12	16	0.0006025796753211821	0.004941153337633693	2.01191375277228e-05	6792.410446962112	6792.410446962112	6792.410446962112	0	0	1	-60	60;
	13	14	0.0007568790811307346	0.006206408465272025	2.527093917830723e-05	6792.410446962112	6792.410446962112	6792.410446962112	0	0	1	-60	60;
	13	5	0.00112366070458674	0.009214017777611268	9.379292554975687e-06	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	13	6	0.001555406933486283	0.012754336854587522	1.2983115465064572e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	14	18	0.001184704521406836	0.009714577075536055	2.224986872573658e-05	5094.307835221583	5094.307835221583	5094.307835221583	0	0	1	-60	60;
	14	3	0.0007624122702652717	0.00625178061617523	4.951870783337223e-05	9473.625097078735	9473.625097078735	9473.625097078735	0	0	1	-60	60;
	14	6	0.00023546671415108084	0.001930827056038863	7.774867939575865e-05	21360.343379262427	21360.343379262427	21360.343379262427	0	0	1	-60	60;
	15	9	0.003025981227480034	0.02481304606533628	2.5258125591245843e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	16	4	0.0005576796429122965	0.004572973071880832	2.4400065523473356e-05	7775.522485338205	7775.522485338205	7775.522485338205	0	0	1	-60	60;
	16	8	0.0022733219892160884	0.018641240311571924	1.8975614188055853e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	17	3	0.0016284779498957418	0.013353519189145081	2.260173289752433e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	17	5	0.0005083825972120511	0.004168737297138818	1.3588646366447507e-05	6077.419873597679	6077.419873597679	6077.419873597679	0	0	1	-60	60;
	18	20	0.00046182876264664046	0.0037869958537024525	3.469431966193496e-05	10188.615670443167	10188.615670443167	10188.615670443167	0	0	1	-60	60;
	19	9	0.00013384008124748897	0.0010974886662294095	7.149917892347253e-05	27169.641787848446	27169.641787848446	27169.641787848446	0	0	1	-60	60;
	20	7	0.00018254458558985453	0.001496865601836807	3.930531396397533e-05	17249.14758241694	17249.14758241694	17249.14758241694	0	0	1	-60	60;
	3	5	0.0014319862400107694	0.01174228716808831	3.8275807873511336e-05	6077.419873597678	6077.419873597678	6077.419873597678	0	0	1	-60	60;
	3	8	0.0011466492608282874	0.009402523938791957	1.591440665318426e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	3	9	0.0005405235524081651	0.004432293129746955	4.8820513613460305e-05	11171.727708819262	11171.727708819262	11171.727708819262	0	0	1	-60	60;
	4	7	0.0018507033179182048	0.01517576720692928	1.5447979786392064e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	4	10	0.0011474237792125823	0.009408874989543172	1.59251562354228e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	7	10	0.00021943617121981208	0.0017993766040024588	3.709096507683684e-05	15282.923505664749	15282.923505664749	15282.923505664749	0	0	1	-60	60;
	8	9	0.0062383406820868985	0.05115439359311257	4.363371135315227e-06	983.1120383760949	983.1120383760949	983.1120383760949	0	0	1	-60	60;
	8	10	0.0005928379554025724	0.0048612712343010945	2.5938341376573478e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
];

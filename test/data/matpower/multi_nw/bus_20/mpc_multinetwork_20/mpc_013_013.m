function mpc = mpc_013_013
%MPC_013_013	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 13 	Weight: 18
%	Time step: 13

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 18;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	2582.986824243746	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5895.170852706913	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	2201.18509186098	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	3701.5099477701083	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	2747.566222137764	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	6754.88820670724	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	394.5912249634915	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3073.5568252619696	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	4838.212846765469	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	5097.670514419961	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	3483.022931669188	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	3170.837744678054	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6074.077071065285	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3739.9472349183243	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	4243.3085611005645	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	5334.399354827476	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	4039.9669499672304	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5902.403764807215	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6386.030517781361	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4664.251171452676	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.875653509283507	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.790453856223932	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.567266613642401	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.57413871918775	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.627033593394879	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.251059168717278	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.704533645204377	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.021534644226614	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.916784853179973	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.823716506259835	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.220475860069623	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.927399019893054	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.60101919348213	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.062080659612338	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.809061347825422	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.0339321896557365	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.469520872739963	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.14973375636686	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.676782987062856	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.814813147131957	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.683463489689265	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.486527869528384	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.058683050747584	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.390855949131488	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.159860383536262	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.172644769520062	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.979927544762287	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.186365044387117	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.143147506110676	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.974977711200212	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.709637136091365	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.673207925451242	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.22179492651718	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.061661265228105	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.591094705741182	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.619760222169097	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.9901909867285448	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.98954649327958	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.557473677202005	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.425105641185931	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.432868807469596	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.377915655703257	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.943501386842723	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.647995712328202	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.755426472378765	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.671832036806968	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	1370.6160143516245	548.2464057406498	3	0.0	268.27481196939215	0.0	0.0	5714.035589507766;	%0 biomass
	2	0.0	0.0	3	0.0	0.4441122419977165	0.0	0.0	6617.754653710759;	%0 offwind
	2	0.0	0.0	3	0.0	0.4337546801127558	0.0	0.0	7502.596735805587;	%0 onwind
	2	0.0	0.0	3	0.0	0.3430594927191341	0.0	0.0	2048.3383455339545;	%0 solar
	2	159892.2502602916	63956.90010411663	3	0.0	815.7767870423038	0.0	0.0	45579.886359587974;	%1 CCGT
	2	195173.9580818642	78069.58323274567	3	0.0	995.7855004176745	0.0	0.0	20756.258012439946;	%1 OCGT
	2	1370.561260976939	548.2245043907756	3	0.0	268.26409492600095	0.0	0.0	5714.0355895077655;	%1 biomass
	2	0.0	0.0	3	0.0	0.44244025185416747	0.0	0.0	7502.596735805587;	%1 onwind
	2	0.0	0.0	3	0.0	0.35891583551289374	0.0	0.0	2048.3383455339545;	%1 solar
	2	140098.98879322596	56039.59551729038	3	0.0	714.7907591491119	0.0	0.0	45579.886359587974;	%10 CCGT
	2	203223.7338939594	81289.49355758376	3	0.0	1036.8557851732621	0.0	0.0	20756.258012439946;	%10 OCGT
	2	1370.5864737841098	548.2345895136439	3	0.0	268.2690299048952	0.0	0.0	5714.0355895077655;	%10 biomass
	2	0.0	0.0	3	0.0	0.4444755937452566	0.0	0.0	6533.240487614452;	%10 offwind
	2	0.0	0.0	3	0.0	0.44639476685151774	0.0	0.0	7502.596735805588;	%10 onwind
	2	0.0	0.0	3	0.0	0.3568010041641963	0.0	0.0	2048.3383455339545;	%10 solar
	2	130726.44225846078	52290.57690338431	3	0.0	666.9716441758203	0.0	0.0	45579.886359587974;	%11 CCGT
	2	194063.48786026845	77625.39514410737	3	0.0	990.1198360217777	0.0	0.0	20756.258012439946;	%11 OCGT
	2	1370.6167855800838	548.2467142320336	3	0.0	268.2749629242677	0.0	0.0	5714.0355895077655;	%11 biomass
	2	0.0	0.0	3	0.0	0.4556559761717418	0.0	0.0	6495.842504400269;	%11 offwind
	2	0.0	0.0	3	0.0	0.4339225853594653	0.0	0.0	7502.596735805587;	%11 onwind
	2	0.0	0.0	3	0.0	0.36171724606846195	0.0	0.0	2048.3383455339545;	%11 solar
	2	193017.66536729448	77207.06614691779	3	0.0	984.7840069759922	0.0	0.0	20756.258012439946;	%12 OCGT
	2	1370.5907361187164	548.2362944474867	3	0.0	268.26986418452077	0.0	0.0	5714.0355895077655;	%12 biomass
	2	0.0	0.0	3	0.0	0.4502871524372878	0.0	0.0	7502.596735805587;	%12 onwind
	2	0.0	0.0	3	0.0	0.3601819303668902	0.0	0.0	2048.3383455339545;	%12 solar
	2	145113.3157597767	58045.32630391068	3	0.0	740.3740599988606	0.0	0.0	45579.886359587974;	%13 CCGT
	2	201859.00094941616	80743.60037976646	3	0.0	1029.892861986817	0.0	0.0	20756.258012439946;	%13 OCGT
	2	1370.5781549310718	548.2312619724287	3	0.0	268.26740163066586	0.0	0.0	5714.0355895077655;	%13 biomass
	2	0.0	0.0	3	0.0	0.4577262974096279	0.0	0.0	7502.596735805587;	%13 onwind
	2	0.0	0.0	3	0.0	0.36438555904202663	0.0	0.0	2048.3383455339545;	%13 solar
	2	149538.5280349143	59815.41121396573	3	0.0	762.9516736475221	0.0	0.0	45579.886359587974;	%14 CCGT
	2	201714.48236205574	80685.7929448223	3	0.0	1029.1555222553866	0.0	0.0	20756.258012439946;	%14 OCGT
	2	1370.5219768963973	548.208790758559	3	0.0	268.25640573427233	0.0	0.0	5714.0355895077655;	%14 biomass
	2	0.0	0.0	3	0.0	0.44504478639807055	0.0	0.0	6841.273529830952;	%14 offwind
	2	0.0	0.0	3	0.0	0.44064004503902804	0.0	0.0	7502.596735805588;	%14 onwind
	2	0.0	0.0	3	0.0	0.35225755520135543	0.0	0.0	2048.3383455339545;	%14 solar
	2	134298.78306822505	53719.513227290015	3	0.0	685.1978727970665	0.0	0.0	45579.886359587974;	%15 CCGT
	2	188366.50834873682	75346.60333949473	3	0.0	961.0536140241675	0.0	0.0	20756.258012439946;	%15 OCGT
	2	1370.5937793396865	548.2375117358746	3	0.0	268.2704598433522	0.0	0.0	5714.0355895077655;	%15 biomass
	2	0.0	0.0	3	0.0	0.45057544961475543	0.0	0.0	7139.55132762114;	%15 offwind
	2	0.0	0.0	3	0.0	0.4460171453031351	0.0	0.0	7502.596735805588;	%15 onwind
	2	0.0	0.0	3	0.0	0.34810962487027586	0.0	0.0	2048.3383455339545;	%15 solar
	2	140580.58634341502	56232.23453736601	3	0.0	717.2478895072194	0.0	0.0	45579.886359587974;	%16 CCGT
	2	194510.4539790089	77804.18159160354	3	0.0	992.4002754031064	0.0	0.0	20756.258012439946;	%16 OCGT
	2	1370.5607525215355	548.2243010086142	3	0.0	268.2639954044892	0.0	0.0	5714.035589507766;	%16 biomass
	2	0.0	0.0	3	0.0	0.4452617374686085	0.0	0.0	7502.596735805587;	%16 onwind
	2	0.0	0.0	3	0.0	0.3580617152902063	0.0	0.0	2048.338345533955;	%16 solar
	2	144535.16808106995	57814.06723242798	3	0.0	737.4243269442343	0.0	0.0	45579.886359587974;	%17 CCGT
	2	189616.4796243489	75846.59184973956	3	0.0	967.4310184915759	0.0	0.0	20756.258012439946;	%17 OCGT
	2	1370.5762299198755	548.2304919679502	3	0.0	268.26702484241054	0.0	0.0	5714.0355895077655;	%17 biomass
	2	0.0	0.0	3	0.0	0.45885805153945314	0.0	0.0	7502.596735805588;	%17 onwind
	2	0.0	0.0	3	0.0	0.3669118463807074	0.0	0.0	2048.338345533955;	%17 solar
	2	151505.43083204384	60602.17233281754	3	0.0	772.9868920002236	0.0	0.0	45579.886359587974;	%18 CCGT
	2	191092.1901611527	76436.87606446109	3	0.0	974.9601538834322	0.0	0.0	20756.258012439946;	%18 OCGT
	2	1370.601037221005	548.240414888402	3	0.0	268.27188045038264	0.0	0.0	5714.0355895077655;	%18 biomass
	2	0.0	0.0	3	0.0	0.4483727151471608	0.0	0.0	7502.596735805587;	%18 onwind
	2	0.0	0.0	3	0.0	0.3599898161454938	0.0	0.0	2048.3383455339545;	%18 solar
	2	141962.70214877123	56785.080859508496	3	0.0	724.2995007590368	0.0	0.0	45579.886359587974;	%19 CCGT
	2	198902.57634884614	79561.03053953845	3	0.0	1014.8090630043168	0.0	0.0	20756.258012439946;	%19 OCGT
	2	1370.569455906085	548.227782362434	3	0.0	268.2656989442327	0.0	0.0	5714.0355895077655;	%19 biomass
	2	0.0	0.0	3	0.0	0.45745264550646075	0.0	0.0	7502.596735805587;	%19 onwind
	2	0.0	0.0	3	0.0	0.35782494885626076	0.0	0.0	2048.3383455339545;	%19 solar
	2	151984.47477330145	60793.78990932058	3	0.0	775.4309937413339	0.0	0.0	45579.886359587974;	%2 CCGT
	2	200866.43654242685	80346.57461697074	3	0.0	1024.8287578695247	0.0	0.0	20756.258012439946;	%2 OCGT
	2	1370.5627530163802	548.225101206552	3	0.0	268.264386967387	0.0	0.0	5714.0355895077655;	%2 biomass
	2	0.0	0.0	3	0.0	0.4508165600704036	0.0	0.0	7502.596735805587;	%2 onwind
	2	0.0	0.0	3	0.0	0.3589743297282216	0.0	0.0	2048.3383455339545;	%2 solar
	2	165742.58799712927	66297.03519885172	3	0.0	845.6254489649453	0.0	0.0	45579.886359587974;	%3 CCGT
	2	206171.4205037455	82468.5682014982	3	0.0	1051.8950025701301	0.0	0.0	20756.258012439946;	%3 OCGT
	2	1370.5227164487267	548.2090865794906	3	0.0	268.25655048908334	0.0	0.0	5714.0355895077655;	%3 biomass
	2	0.0	0.0	3	0.0	0.44593752776835877	0.0	0.0	7107.763863575269;	%3 offwind
	2	0.0	0.0	3	0.0	0.45876450713866357	0.0	0.0	7502.596735805588;	%3 onwind
	2	0.0	0.0	3	0.0	0.3594627431008836	0.0	0.0	2048.3383455339545;	%3 solar
	2	144479.97545914195	57791.99018365678	3	0.0	737.1427319343977	0.0	0.0	45579.886359587974;	%4 CCGT
	2	197588.69733365162	79035.47893346066	3	0.0	1008.1055986410797	0.0	0.0	20756.258012439946;	%4 OCGT
	2	1370.6041202311615	548.2416480924645	3	0.0	268.27248389727174	0.0	0.0	5714.0355895077655;	%4 biomass
	2	0.0	0.0	3	0.0	0.4487365372808509	0.0	0.0	7502.596735805587;	%4 onwind
	2	0.0	0.0	3	0.0	0.3518016957840727	0.0	0.0	2048.3383455339545;	%4 solar
	2	148399.7647726997	59359.90590907988	3	0.0	757.1416570035699	0.0	0.0	45579.886359587974;	%5 CCGT
	2	194028.3585530218	77611.34342120872	3	0.0	989.9406048623559	0.0	0.0	20756.258012439946;	%5 OCGT
	2	1370.5858845045889	548.2343538018356	3	0.0	268.2689145634349	0.0	0.0	5714.0355895077655;	%5 biomass
	2	0.0	0.0	3	0.0	0.4432810928034283	0.0	0.0	7502.596735805587;	%5 onwind
	2	0.0	0.0	3	0.0	0.3603060982425506	0.0	0.0	2048.3383455339545;	%5 solar
	2	145504.5127412539	58201.80509650156	3	0.0	742.3699629655811	0.0	0.0	45579.886359587974;	%6 CCGT
	2	200958.51275364406	80383.40510145763	3	0.0	1025.2985344573676	0.0	0.0	20756.258012439946;	%6 OCGT
	2	1370.608611581753	548.2434446327012	3	0.0	268.27336300288766	0.0	0.0	5714.0355895077655;	%6 biomass
	2	0.0	0.0	3	0.0	0.4501744435043924	0.0	0.0	7502.596735805587;	%6 onwind
	2	0.0	0.0	3	0.0	0.3675394242711416	0.0	0.0	2048.3383455339545;	%6 solar
	2	137917.2129277718	55166.88517110871	3	0.0	703.6592496314887	0.0	0.0	45579.886359587974;	%7 CCGT
	2	200099.37317877295	80039.74927150918	3	0.0	1020.9151692794537	0.0	0.0	20756.258012439946;	%7 OCGT
	2	1370.607143845425	548.24285753817	3	0.0	268.2730757184234	0.0	0.0	5714.0355895077655;	%7 biomass
	2	0.0	0.0	3	0.0	0.45509841967458814	0.0	0.0	7502.596735805587;	%7 onwind
	2	0.0	0.0	3	0.0	0.35910187357913936	0.0	0.0	2048.338345533955;	%7 solar
	2	155934.90170984276	62373.9606839371	3	0.0	795.5862332134833	0.0	0.0	45579.886359587974;	%8 CCGT
	2	204489.97689759324	81795.9907590373	3	0.0	1043.31620866119	0.0	0.0	20756.258012439946;	%8 OCGT
	2	1370.6097492823553	548.243899712942	3	0.0	268.27358568846256	0.0	0.0	5714.0355895077655;	%8 biomass
	2	0.0	0.0	3	0.0	0.4600381453745848	0.0	0.0	7502.596735805587;	%8 onwind
	2	0.0	0.0	3	0.0	0.35864349991615874	0.0	0.0	2048.3383455339545;	%8 solar
	2	149528.35009812808	59811.340039251234	3	0.0	762.8997453986126	0.0	0.0	45579.886359587974;	%9 CCGT
	2	201754.91538819208	80701.96615527684	3	0.0	1029.3618132050615	0.0	0.0	20756.258012439946;	%9 OCGT
	2	1370.5559591557394	548.2223836622958	3	0.0	268.2630571845252	0.0	0.0	5714.0355895077655;	%9 biomass
	2	0.0	0.0	3	0.0	0.4474031792164636	0.0	0.0	7502.596735805587;	%9 onwind
	2	0.0	0.0	3	0.0	0.3628234052683575	0.0	0.0	2048.3383455339545;	%9 solar
	2	0.0	0.0	3	0.0	0.17402556467785107	0.0	0.0	364.40797935510903;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.13051917350838832	0.0	0.0	364.40797935510903;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1895975350097391	0.0	0.0	364.40797935510903;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.14219815125730434	0.0	0.0	364.40797935510903;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19783486091573668	0.0	0.0	364.40797935510903;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.14837614568680252	0.0	0.0	364.40797935510903;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16521199371955264	0.0	0.0	364.40797935510903;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.12390899528966448	0.0	0.0	364.40797935510903;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19655436058402223	0.0	0.0	364.40797935510903;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.14741577043801668	0.0	0.0	364.40797935510903;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.189085427922496	0.0	0.0	364.40797935510903;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.14181407094187198	0.0	0.0	364.40797935510903;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17998422808830444	0.0	0.0	364.40797935510903;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.13498817106622835	0.0	0.0	364.40797935510903;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18269747967870908	0.0	0.0	364.40797935510903;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.1806411585271795	0.0	0.0	364.40797935510903;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.13548086889538463	0.0	0.0	364.40797935510903;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17590383582049354	0.0	0.0	364.40797935510903;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.13192787686537016	0.0	0.0	364.40797935510903;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1726706667576203	0.0	0.0	364.40797935510903;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.12950300006821522	0.0	0.0	364.40797935510903;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16535745522795475	0.0	0.0	364.40797935510903;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.12401809142096606	0.0	0.0	364.40797935510903;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1823841375153478	0.0	0.0	364.40797935510903;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.13678810313651085	0.0	0.0	364.40797935510903;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16741972662667232	0.0	0.0	364.40797935510903;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.12556479497000422	0.0	0.0	364.40797935510903;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16842695804734656	0.0	0.0	364.40797935510903;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.12632021853550993	0.0	0.0	364.40797935510903;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1951854315916587	0.0	0.0	364.40797935510903;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.14638907369374402	0.0	0.0	364.40797935510903;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.17867711344124457	0.0	0.0	364.40797935510903;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.13400783508093342	0.0	0.0	364.40797935510903;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1857004664754963	0.0	0.0	364.40797935510903;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.1392753498566222	0.0	0.0	364.40797935510903;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18826276094834357	0.0	0.0	364.40797935510903;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.14119707071125767	0.0	0.0	364.40797935510903;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16715507956350417	0.0	0.0	364.40797935510903;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.19456673829673624	0.0	0.0	364.40797935510903;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.14592505372255218	0.0	0.0	364.40797935510903;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_010_013
%MPC_010_013	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 10 	Weight: 12
%	Time step: 13

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 12;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1247.4633870453981	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4929.267870850888	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	621.5783790578429	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-1178.0676576009087	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	-266.09549790851923	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	6102.291667593689	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	317.65152484843554	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2599.482365752522	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2293.648349522027	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	5350.811107266392	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	-1363.216557052007	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	4572.339920177643	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4865.19569123897	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4721.349515833565	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2262.5954369501533	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	2877.885007952312	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	7294.334545680762	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4930.146179398554	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3837.1957370511086	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3810.325265639037	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.454043106386365	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.458850065145747	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.175311335990926	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.4103238679908183	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.256960180508212	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.637553099036712	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.930419572192546	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.710435782440118	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.240340794987949	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.95209426875109	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.451684630155262	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.22062894656568846	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.157966491462918	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.35245973786503904	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.791085996248574	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.21611898468116	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.56443519189997	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.3876428100106	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.86628265457347	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.614542470073138	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.823725230608005	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.4583526721869973	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.320637160007635	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.403723730125547	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.628865827151229	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.236632217917439	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.821008478404124	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.576791998863186	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5848635478670547	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.459060886228825	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.953426044447145	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.242456277566584	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.387465154067363	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.2314395105956937	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.238094168136675	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.658923816681389	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.023966874271881	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.373626601781686	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.957586037188482	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.54359586922081	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.748582180373994	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.50936209372626	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.748560845415101	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.9069299847414865	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	913.7440095677497	365.49760382709985	3	0.0	178.84987464626144	0.0	0.0	3809.3570596718446;	%0 biomass
	2	0.0	0.0	3	0.0	0.2960748279984777	0.0	0.0	4411.836435807173;	%0 offwind
	2	0.0	0.0	3	0.0	0.28916978674183724	0.0	0.0	5001.731157203724;	%0 onwind
	2	0.0	0.0	3	0.0	0.22870632847942274	0.0	0.0	1365.5588970226363;	%0 solar
	2	106594.83350686106	42637.93340274442	3	0.0	543.8511913615359	0.0	0.0	30386.590906391983;	%1 CCGT
	2	130115.97205457614	52046.38882183045	3	0.0	663.8570002784496	0.0	0.0	13837.505341626631;	%1 OCGT
	2	913.7075073179592	365.48300292718375	3	0.0	178.84272995066732	0.0	0.0	3809.3570596718437;	%1 biomass
	2	0.0	0.0	3	0.0	0.29496016790277835	0.0	0.0	5001.731157203724;	%1 onwind
	2	0.0	0.0	3	0.0	0.2392772236752625	0.0	0.0	1365.5588970226363;	%1 solar
	2	93399.32586215064	37359.73034486025	3	0.0	476.5271727660746	0.0	0.0	30386.590906391983;	%10 CCGT
	2	135482.4892626396	54192.995705055844	3	0.0	691.2371901155082	0.0	0.0	13837.505341626631;	%10 OCGT
	2	913.7243158560732	365.4897263424292	3	0.0	178.8460199365968	0.0	0.0	3809.3570596718437;	%10 biomass
	2	0.0	0.0	3	0.0	0.2963170624968377	0.0	0.0	4355.493658409635;	%10 offwind
	2	0.0	0.0	3	0.0	0.29759651123434516	0.0	0.0	5001.731157203725;	%10 onwind
	2	0.0	0.0	3	0.0	0.2378673361094642	0.0	0.0	1365.5588970226363;	%10 solar
	2	87150.96150564053	34860.38460225621	3	0.0	444.64776278388024	0.0	0.0	30386.590906391983;	%11 CCGT
	2	129375.65857351231	51750.26342940492	3	0.0	660.0798906811851	0.0	0.0	13837.505341626631;	%11 OCGT
	2	913.7445237200559	365.49780948802237	3	0.0	178.84997528284515	0.0	0.0	3809.3570596718437;	%11 biomass
	2	0.0	0.0	3	0.0	0.30377065078116117	0.0	0.0	4330.56166960018;	%11 offwind
	2	0.0	0.0	3	0.0	0.2892817235729769	0.0	0.0	5001.731157203724;	%11 onwind
	2	0.0	0.0	3	0.0	0.24114483071230797	0.0	0.0	1365.5588970226363;	%11 solar
	2	128678.44357819633	51471.377431278524	3	0.0	656.5226713173281	0.0	0.0	13837.505341626631;	%12 OCGT
	2	913.7271574124777	365.4908629649911	3	0.0	178.84657612301385	0.0	0.0	3809.3570596718437;	%12 biomass
	2	0.0	0.0	3	0.0	0.3001914349581919	0.0	0.0	5001.731157203724;	%12 onwind
	2	0.0	0.0	3	0.0	0.24012128691126014	0.0	0.0	1365.5588970226363;	%12 solar
	2	96742.2105065178	38696.88420260712	3	0.0	493.5827066659071	0.0	0.0	30386.590906391983;	%13 CCGT
	2	134572.6672996108	53829.06691984431	3	0.0	686.5952413245448	0.0	0.0	13837.505341626631;	%13 OCGT
	2	913.718769954048	365.48750798161916	3	0.0	178.8449344204439	0.0	0.0	3809.3570596718437;	%13 biomass
	2	0.0	0.0	3	0.0	0.30515086493975196	0.0	0.0	5001.731157203724;	%13 onwind
	2	0.0	0.0	3	0.0	0.24292370602801774	0.0	0.0	1365.5588970226363;	%13 solar
	2	99692.35202327621	39876.94080931049	3	0.0	508.634449098348	0.0	0.0	30386.590906391983;	%14 CCGT
	2	134476.32157470385	53790.52862988153	3	0.0	686.103681503591	0.0	0.0	13837.505341626631;	%14 OCGT
	2	913.6813179309315	365.4725271723726	3	0.0	178.83760382284822	0.0	0.0	3809.3570596718437;	%14 biomass
	2	0.0	0.0	3	0.0	0.29669652426538035	0.0	0.0	4560.8490198873005;	%14 offwind
	2	0.0	0.0	3	0.0	0.2937600300260187	0.0	0.0	5001.731157203725;	%14 onwind
	2	0.0	0.0	3	0.0	0.23483837013423695	0.0	0.0	1365.5588970226363;	%14 solar
	2	89532.52204548336	35813.008818193346	3	0.0	456.79858186471097	0.0	0.0	30386.590906391983;	%15 CCGT
	2	125577.67223249121	50231.06889299649	3	0.0	640.702409349445	0.0	0.0	13837.505341626631;	%15 OCGT
	2	913.7291862264576	365.4916744905831	3	0.0	178.8469732289015	0.0	0.0	3809.3570596718437;	%15 biomass
	2	0.0	0.0	3	0.0	0.3003836330765036	0.0	0.0	4759.700885080761;	%15 offwind
	2	0.0	0.0	3	0.0	0.2973447635354234	0.0	0.0	5001.731157203725;	%15 onwind
	2	0.0	0.0	3	0.0	0.23207308324685055	0.0	0.0	1365.5588970226363;	%15 solar
	2	93720.39089561	37488.156358244	3	0.0	478.1652596714796	0.0	0.0	30386.590906391983;	%16 CCGT
	2	129673.63598600592	51869.45439440236	3	0.0	661.600183602071	0.0	0.0	13837.505341626631;	%16 OCGT
	2	913.7071683476904	365.48286733907617	3	0.0	178.84266360299281	0.0	0.0	3809.3570596718446;	%16 biomass
	2	0.0	0.0	3	0.0	0.29684115831240565	0.0	0.0	5001.731157203724;	%16 onwind
	2	0.0	0.0	3	0.0	0.23870781019347087	0.0	0.0	1365.5588970226365;	%16 solar
	2	96356.7787207133	38542.71148828532	3	0.0	491.6162179628229	0.0	0.0	30386.590906391983;	%17 CCGT
	2	126410.9864162326	50564.394566493036	3	0.0	644.9540123277172	0.0	0.0	13837.505341626631;	%17 OCGT
	2	913.7174866132503	365.4869946453001	3	0.0	178.8446832282737	0.0	0.0	3809.3570596718437;	%17 biomass
	2	0.0	0.0	3	0.0	0.3059053676929688	0.0	0.0	5001.731157203725;	%17 onwind
	2	0.0	0.0	3	0.0	0.24460789758713825	0.0	0.0	1365.5588970226365;	%17 solar
	2	101003.6205546959	40401.44822187836	3	0.0	515.3245946668158	0.0	0.0	30386.590906391983;	%18 CCGT
	2	127394.79344076848	50957.91737630739	3	0.0	649.9734359222881	0.0	0.0	13837.505341626631;	%18 OCGT
	2	913.7340248140033	365.4936099256014	3	0.0	178.8479203002551	0.0	0.0	3809.3570596718437;	%18 biomass
	2	0.0	0.0	3	0.0	0.29891514343144054	0.0	0.0	5001.731157203724;	%18 onwind
	2	0.0	0.0	3	0.0	0.23999321076366253	0.0	0.0	1365.5588970226363;	%18 solar
	2	94641.80143251416	37856.72057300567	3	0.0	482.86633383935794	0.0	0.0	30386.590906391983;	%19 CCGT
	2	132601.71756589742	53040.687026358966	3	0.0	676.5393753362113	0.0	0.0	13837.505341626631;	%19 OCGT
	2	913.7129706040566	365.4851882416227	3	0.0	178.84379929615517	0.0	0.0	3809.3570596718437;	%19 biomass
	2	0.0	0.0	3	0.0	0.30496843033764054	0.0	0.0	5001.731157203724;	%19 onwind
	2	0.0	0.0	3	0.0	0.23854996590417382	0.0	0.0	1365.5588970226363;	%19 solar
	2	101322.98318220096	40529.19327288039	3	0.0	516.9539958275559	0.0	0.0	30386.590906391983;	%2 CCGT
	2	133910.95769495124	53564.3830779805	3	0.0	683.2191719130165	0.0	0.0	13837.505341626631;	%2 OCGT
	2	913.70850201092	365.48340080436805	3	0.0	178.84292464492466	0.0	0.0	3809.3570596718437;	%2 biomass
	2	0.0	0.0	3	0.0	0.30054437338026907	0.0	0.0	5001.731157203724;	%2 onwind
	2	0.0	0.0	3	0.0	0.2393162198188144	0.0	0.0	1365.5588970226363;	%2 solar
	2	110495.05866475285	44198.02346590115	3	0.0	563.7502993099636	0.0	0.0	30386.590906391983;	%3 CCGT
	2	137447.61366916366	54979.04546766546	3	0.0	701.2633350467534	0.0	0.0	13837.505341626631;	%3 OCGT
	2	913.6818109658177	365.4727243863271	3	0.0	178.83770032605554	0.0	0.0	3809.3570596718437;	%3 biomass
	2	0.0	0.0	3	0.0	0.2972916851789058	0.0	0.0	4738.509242383513;	%3 offwind
	2	0.0	0.0	3	0.0	0.305843004759109	0.0	0.0	5001.731157203725;	%3 onwind
	2	0.0	0.0	3	0.0	0.2396418287339224	0.0	0.0	1365.5588970226363;	%3 solar
	2	96319.98363942797	38527.99345577119	3	0.0	491.42848795626514	0.0	0.0	30386.590906391983;	%4 CCGT
	2	131725.79822243442	52690.31928897377	3	0.0	672.0703990940532	0.0	0.0	13837.505341626631;	%4 OCGT
	2	913.7360801541076	365.494432061643	3	0.0	178.84832259818117	0.0	0.0	3809.3570596718437;	%4 biomass
	2	0.0	0.0	3	0.0	0.29915769152056726	0.0	0.0	5001.731157203724;	%4 onwind
	2	0.0	0.0	3	0.0	0.23453446385604848	0.0	0.0	1365.5588970226363;	%4 solar
	2	98933.17651513313	39573.27060605325	3	0.0	504.7611046690466	0.0	0.0	30386.590906391983;	%5 CCGT
	2	129352.23903534785	51740.89561413914	3	0.0	659.9604032415707	0.0	0.0	13837.505341626631;	%5 OCGT
	2	913.7239230030593	365.48956920122373	3	0.0	178.84594304228995	0.0	0.0	3809.3570596718437;	%5 biomass
	2	0.0	0.0	3	0.0	0.2955207285356189	0.0	0.0	5001.731157203724;	%5 onwind
	2	0.0	0.0	3	0.0	0.24020406549503373	0.0	0.0	1365.5588970226363;	%5 solar
	2	97003.00849416925	38801.203397667705	3	0.0	494.9133086437207	0.0	0.0	30386.590906391983;	%6 CCGT
	2	133972.34183576272	53588.93673430508	3	0.0	683.5323563049118	0.0	0.0	13837.505341626631;	%6 OCGT
	2	913.7390743878353	365.49562975513413	3	0.0	178.84890866859178	0.0	0.0	3809.3570596718437;	%6 biomass
	2	0.0	0.0	3	0.0	0.3001162956695949	0.0	0.0	5001.731157203724;	%6 onwind
	2	0.0	0.0	3	0.0	0.24502628284742772	0.0	0.0	1365.5588970226363;	%6 solar
	2	91944.80861851452	36777.92344740581	3	0.0	469.10616642099245	0.0	0.0	30386.590906391983;	%7 CCGT
	2	133399.58211918199	53359.832847672784	3	0.0	680.6101128529691	0.0	0.0	13837.505341626631;	%7 OCGT
	2	913.7380958969501	365.49523835878006	3	0.0	178.84871714561558	0.0	0.0	3809.3570596718437;	%7 biomass
	2	0.0	0.0	3	0.0	0.3033989464497254	0.0	0.0	5001.731157203724;	%7 onwind
	2	0.0	0.0	3	0.0	0.23940124905275958	0.0	0.0	1365.5588970226365;	%7 solar
	2	103956.60113989517	41582.64045595806	3	0.0	530.3908221423221	0.0	0.0	30386.590906391983;	%8 CCGT
	2	136326.65126506216	54530.66050602487	3	0.0	695.5441391074601	0.0	0.0	13837.505341626631;	%8 OCGT
	2	913.7398328549035	365.49593314196136	3	0.0	178.8490571256417	0.0	0.0	3809.3570596718437;	%8 biomass
	2	0.0	0.0	3	0.0	0.3066920969163899	0.0	0.0	5001.731157203724;	%8 onwind
	2	0.0	0.0	3	0.0	0.23909566661077247	0.0	0.0	1365.5588970226363;	%8 solar
	2	99685.5667320854	39874.226692834156	3	0.0	508.59983026574173	0.0	0.0	30386.590906391983;	%9 CCGT
	2	134503.27692546137	53801.31077018456	3	0.0	686.2412088033743	0.0	0.0	13837.505341626631;	%9 OCGT
	2	913.7039727704929	365.48158910819717	3	0.0	178.84203812301683	0.0	0.0	3809.3570596718437;	%9 biomass
	2	0.0	0.0	3	0.0	0.2982687861443091	0.0	0.0	5001.731157203724;	%9 onwind
	2	0.0	0.0	3	0.0	0.24188227017890498	0.0	0.0	1365.5588970226363;	%9 solar
	2	0.0	0.0	3	0.0	0.11601704311856739	0.0	0.0	242.93865290340602;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.08701278233892554	0.0	0.0	242.93865290340602;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1263983566731594	0.0	0.0	242.93865290340602;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.09479876750486955	0.0	0.0	242.93865290340602;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13188990727715777	0.0	0.0	242.93865290340602;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.09891743045786833	0.0	0.0	242.93865290340602;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11014132914636843	0.0	0.0	242.93865290340602;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.08260599685977632	0.0	0.0	242.93865290340602;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13103624038934814	0.0	0.0	242.93865290340602;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.0982771802920111	0.0	0.0	242.93865290340602;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12605695194833066	0.0	0.0	242.93865290340602;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.094542713961248	0.0	0.0	242.93865290340602;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11998948539220296	0.0	0.0	242.93865290340602;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.08999211404415222	0.0	0.0	242.93865290340602;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12179831978580605	0.0	0.0	242.93865290340602;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.12042743901811966	0.0	0.0	242.93865290340602;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.09032057926358975	0.0	0.0	242.93865290340602;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11726922388032904	0.0	0.0	242.93865290340602;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.08795191791024679	0.0	0.0	242.93865290340602;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11511377783841353	0.0	0.0	242.93865290340602;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.08633533337881015	0.0	0.0	242.93865290340602;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11023830348530317	0.0	0.0	242.93865290340602;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.08267872761397738	0.0	0.0	242.93865290340602;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12158942501023186	0.0	0.0	242.93865290340602;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.0911920687576739	0.0	0.0	242.93865290340602;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11161315108444822	0.0	0.0	242.93865290340602;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.08370986331333616	0.0	0.0	242.93865290340602;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11228463869823105	0.0	0.0	242.93865290340602;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.0842134790236733	0.0	0.0	242.93865290340602;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1301236210611058	0.0	0.0	242.93865290340602;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.09759271579582936	0.0	0.0	242.93865290340602;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11911807562749638	0.0	0.0	242.93865290340602;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.08933855672062228	0.0	0.0	242.93865290340602;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1238003109836642	0.0	0.0	242.93865290340602;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.09285023323774816	0.0	0.0	242.93865290340602;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12550850729889573	0.0	0.0	242.93865290340602;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.0941313804741718	0.0	0.0	242.93865290340602;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11143671970900278	0.0	0.0	242.93865290340602;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.12971115886449083	0.0	0.0	242.93865290340602;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.09728336914836813	0.0	0.0	242.93865290340602;	%DE1 76 PHS (pump mode)
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

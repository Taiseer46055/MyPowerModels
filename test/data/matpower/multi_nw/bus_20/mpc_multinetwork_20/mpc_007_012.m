function mpc = mpc_007_012
%MPC_007_012	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 7 	Weight: 51
%	Time step: 12

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 51;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	-212.9645594059422	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6257.984836294663	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1950.4679150684242	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	181.52975064083913	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	599.4823723840166	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	7632.272347849627	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	334.592383042918	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3255.299913737562	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	1927.1820299462486	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	4365.767335842915	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1851.0367514893394	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	6769.16146577797	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6442.085865902664	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4214.126047125584	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2658.320496543161	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	443.4042457837128	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	9123.833045727055	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5587.939528953177	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3519.3863115731124	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4753.415289911124	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6239020904240467	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.182551975698907	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.38549413541216937	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.010786279237273875	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7663724761847086	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.38295457302607966	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.19101704990419507	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.4738984675066633	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6230435302065913	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.417027344577657	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.13654496574851635	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.53250027123958	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.1269261566268631	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.426633113578477	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5077870378383472	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.023197935126846064	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8830701284548725	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.4597745855638876	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.769418759422392	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.1166131698425599	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.09409552266985	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.7190099767539486	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.178721324605671	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.073471761216362	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.514805080151184	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.744468177903147	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.34908489099318524	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.661042556443018	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.2426574836741189	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.508414123454296	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.050067947491147	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.279674686851503	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.007854639658398428	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.685198467192375	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.1876364711730325	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	3883.4120406629363	1553.3648162651743	3	0.0	760.1119672466111	0.0	0.0	16189.767503605339;	%0 biomass
	2	0.0	0.0	3	0.0	1.2583180189935301	0.0	0.0	18750.304852180485;	%0 offwind
	2	0.0	0.0	3	0.0	1.228971593652808	0.0	0.0	21257.35741811583;	%0 onwind
	2	0.0	0.0	3	0.0	0.9720018960375467	0.0	0.0	5803.625312346205;	%0 solar
	2	453028.0424041595	181211.21696166377	3	0.0	2311.3675632865275	0.0	0.0	129143.01135216594;	%1 CCGT
	2	552992.8812319486	221197.15249277942	3	0.0	2821.3922511834107	0.0	0.0	58809.397701913185;	%1 OCGT
	2	3883.256906101327	1553.3027624405308	3	0.0	760.081602290336	0.0	0.0	16189.767503605335;	%1 biomass
	2	0.0	0.0	3	0.0	1.2535807135868078	0.0	0.0	21257.35741811583;	%1 onwind
	2	0.0	0.0	3	0.0	1.0169282006198657	0.0	0.0	5803.625312346205;	%1 solar
	2	396947.13491414016	158778.85396565608	3	0.0	2025.240484255817	0.0	0.0	129143.01135216594;	%10 CCGT
	2	575800.5793662183	230320.23174648732	3	0.0	2937.7580579909095	0.0	0.0	58809.397701913185;	%10 OCGT
	2	3883.328342388311	1553.3313369553243	3	0.0	760.0955847305364	0.0	0.0	16189.767503605335;	%10 biomass
	2	0.0	0.0	3	0.0	1.2593475156115603	0.0	0.0	18510.848048240947;	%10 offwind
	2	0.0	0.0	3	0.0	1.264785172745967	0.0	0.0	21257.357418115833;	%10 onwind
	2	0.0	0.0	3	0.0	1.0109361784652229	0.0	0.0	5803.625312346205;	%10 solar
	2	370391.5863989722	148156.63455958889	3	0.0	1889.7529918314908	0.0	0.0	129143.01135216594;	%11 CCGT
	2	549846.5489374272	219938.61957497091	3	0.0	2805.339535395037	0.0	0.0	58809.397701913185;	%11 OCGT
	2	3883.414225810238	1553.365690324095	3	0.0	760.1123949520919	0.0	0.0	16189.767503605335;	%11 biomass
	2	0.0	0.0	3	0.0	1.291025265819935	0.0	0.0	18404.88709580076;	%11 offwind
	2	0.0	0.0	3	0.0	1.2294473251851517	0.0	0.0	21257.35741811583;	%11 onwind
	2	0.0	0.0	3	0.0	1.0248655305273089	0.0	0.0	5803.625312346205;	%11 solar
	2	546883.3852073344	218753.35408293374	3	0.0	2790.221353098644	0.0	0.0	58809.397701913185;	%12 OCGT
	2	3883.34041900303	1553.3361676012123	3	0.0	760.0979485228088	0.0	0.0	16189.767503605335;	%12 biomass
	2	0.0	0.0	3	0.0	1.2758135985723154	0.0	0.0	21257.35741811583;	%12 onwind
	2	0.0	0.0	3	0.0	1.0205154693728555	0.0	0.0	5803.625312346205;	%12 solar
	2	411154.3946527007	164461.75786108026	3	0.0	2097.726503330105	0.0	0.0	129143.01135216594;	%13 CCGT
	2	571933.8360233458	228773.53440933832	3	0.0	2918.0297756293153	0.0	0.0	58809.397701913185;	%13 OCGT
	2	3883.3047723047034	1553.3219089218815	3	0.0	760.0909712868865	0.0	0.0	16189.767503605335;	%13 biomass
	2	0.0	0.0	3	0.0	1.2968911759939457	0.0	0.0	21257.35741811583;	%13 onwind
	2	0.0	0.0	3	0.0	1.0324257506190755	0.0	0.0	5803.625312346205;	%13 solar
	2	423692.4960989239	169476.99843956958	3	0.0	2161.696408667979	0.0	0.0	129143.01135216594;	%14 CCGT
	2	571524.3666924913	228609.74667699652	3	0.0	2915.940646390262	0.0	0.0	58809.397701913185;	%14 OCGT
	2	3883.145601206459	1553.2582404825837	3	0.0	760.0598162471049	0.0	0.0	16189.767503605335;	%14 biomass
	2	0.0	0.0	3	0.0	1.2609602281278665	0.0	0.0	19383.608334521028;	%14 offwind
	2	0.0	0.0	3	0.0	1.2484801276105795	0.0	0.0	21257.357418115833;	%14 onwind
	2	0.0	0.0	3	0.0	0.998063073070507	0.0	0.0	5803.625312346205;	%14 solar
	2	380513.2186933043	152205.28747732172	3	0.0	1941.3939729250217	0.0	0.0	129143.01135216594;	%15 CCGT
	2	533705.1069880876	213482.0427952351	3	0.0	2722.9852397351415	0.0	0.0	58809.397701913185;	%15 OCGT
	2	3883.3490414624453	1553.339616584978	3	0.0	760.0996362228314	0.0	0.0	16189.767503605335;	%15 biomass
	2	0.0	0.0	3	0.0	1.2766304405751403	0.0	0.0	20228.728761593233;	%15 offwind
	2	0.0	0.0	3	0.0	1.2637152450255495	0.0	0.0	21257.357418115833;	%15 onwind
	2	0.0	0.0	3	0.0	0.9863106037991148	0.0	0.0	5803.625312346205;	%15 solar
	2	398311.6613063426	159324.664522537	3	0.0	2032.2023536037882	0.0	0.0	129143.01135216594;	%16 CCGT
	2	551112.9529405251	220445.18117621006	3	0.0	2811.800780308801	0.0	0.0	58809.397701913185;	%16 OCGT
	2	3883.255465477684	1553.3021861910736	3	0.0	760.0813203127194	0.0	0.0	16189.767503605339;	%16 biomass
	2	0.0	0.0	3	0.0	1.261574922827724	0.0	0.0	21257.35741811583;	%16 onwind
	2	0.0	0.0	3	0.0	1.014508193322251	0.0	0.0	5803.625312346206;	%16 solar
	2	409516.3095630315	163806.52382521261	3	0.0	2089.3689263419974	0.0	0.0	129143.01135216594;	%17 CCGT
	2	537246.6922689886	214898.67690759542	3	0.0	2741.0545523927985	0.0	0.0	58809.397701913185;	%17 OCGT
	2	3883.299318106314	1553.3197272425255	3	0.0	760.0899037201632	0.0	0.0	16189.767503605335;	%17 biomass
	2	0.0	0.0	3	0.0	1.3000978126951173	0.0	0.0	21257.357418115833;	%17 onwind
	2	0.0	0.0	3	0.0	1.0395835647453375	0.0	0.0	5803.625312346206;	%17 solar
	2	429265.38735745754	171706.15494298303	3	0.0	2190.129527333967	0.0	0.0	129143.01135216594;	%18 CCGT
	2	541427.872123266	216571.14884930642	3	0.0	2762.3871026697243	0.0	0.0	58809.397701913185;	%18 OCGT
	2	3883.3696054595143	1553.3478421838058	3	0.0	760.1036612760842	0.0	0.0	16189.767503605335;	%18 biomass
	2	0.0	0.0	3	0.0	1.2703893595836222	0.0	0.0	21257.35741811583;	%18 onwind
	2	0.0	0.0	3	0.0	1.0199711457455658	0.0	0.0	5803.625312346205;	%18 solar
	2	402227.6560881852	160891.06243527407	3	0.0	2052.1819188172713	0.0	0.0	129143.01135216594;	%19 CCGT
	2	563557.299655064	225422.9198620256	3	0.0	2875.292345178898	0.0	0.0	58809.397701913185;	%19 OCGT
	2	3883.280125067241	1553.3120500268965	3	0.0	760.0861470086594	0.0	0.0	16189.767503605335;	%19 biomass
	2	0.0	0.0	3	0.0	1.2961158289349721	0.0	0.0	21257.35741811583;	%19 onwind
	2	0.0	0.0	3	0.0	1.0138373550927389	0.0	0.0	5803.625312346205;	%19 solar
	2	430622.67852435407	172249.07140974165	3	0.0	2197.0544822671127	0.0	0.0	129143.01135216594;	%2 CCGT
	2	569121.5702035427	227648.62808141712	3	0.0	2903.6814806303196	0.0	0.0	58809.397701913185;	%2 OCGT
	2	3883.2611335464103	1553.3044534185642	3	0.0	760.0824297409299	0.0	0.0	16189.767503605335;	%2 biomass
	2	0.0	0.0	3	0.0	1.2773135868661436	0.0	0.0	21257.35741811583;	%2 onwind
	2	0.0	0.0	3	0.0	1.0170939342299612	0.0	0.0	5803.625312346205;	%2 solar
	2	469603.99932519963	187841.59973007985	3	0.0	2395.938772067345	0.0	0.0	129143.01135216594;	%3 CCGT
	2	584152.3580939457	233660.94323757823	3	0.0	2980.3691739487017	0.0	0.0	58809.397701913185;	%3 OCGT
	2	3883.1476966047253	1553.2590786418903	3	0.0	760.0602263857361	0.0	0.0	16189.767503605335;	%3 biomass
	2	0.0	0.0	3	0.0	1.2634896620103497	0.0	0.0	20138.66428012993;	%3 offwind
	2	0.0	0.0	3	0.0	1.2998327702262134	0.0	0.0	21257.357418115833;	%3 onwind
	2	0.0	0.0	3	0.0	1.0184777721191702	0.0	0.0	5803.625312346205;	%3 solar
	2	409359.93046756886	163743.97218702754	3	0.0	2088.5710738141265	0.0	0.0	129143.01135216594;	%4 CCGT
	2	559834.6424453463	223933.85697813853	3	0.0	2856.299196149726	0.0	0.0	58809.397701913185;	%4 OCGT
	2	3883.3783406549574	1553.3513362619829	3	0.0	760.10537104227	0.0	0.0	16189.767503605335;	%4 biomass
	2	0.0	0.0	3	0.0	1.2714201889624108	0.0	0.0	21257.35741811583;	%4 onwind
	2	0.0	0.0	3	0.0	0.9967714713882061	0.0	0.0	5803.625312346205;	%4 solar
	2	420466.0001893158	168186.40007572633	3	0.0	2145.234694843448	0.0	0.0	129143.01135216594;	%5 CCGT
	2	549747.0159002284	219898.80636009134	3	0.0	2804.8317137766753	0.0	0.0	58809.397701913185;	%5 OCGT
	2	3883.326672763002	1553.3306691052007	3	0.0	760.0952579297323	0.0	0.0	16189.767503605335;	%5 biomass
	2	0.0	0.0	3	0.0	1.2559630962763804	0.0	0.0	21257.35741811583;	%5 onwind
	2	0.0	0.0	3	0.0	1.0208672783538935	0.0	0.0	5803.625312346205;	%5 solar
	2	412262.78610021935	164905.11444008772	3	0.0	2103.3815617358127	0.0	0.0	129143.01135216594;	%6 CCGT
	2	569382.4528019915	227752.9811207966	3	0.0	2905.012514295875	0.0	0.0	58809.397701913185;	%6 OCGT
	2	3883.3910661483	1553.3564264593201	3	0.0	760.1078618415149	0.0	0.0	16189.767503605335;	%6 biomass
	2	0.0	0.0	3	0.0	1.2754942565957785	0.0	0.0	21257.35741811583;	%6 onwind
	2	0.0	0.0	3	0.0	1.041361702101568	0.0	0.0	5803.625312346205;	%6 solar
	2	390765.43662868673	156306.1746514747	3	0.0	1993.701207289218	0.0	0.0	129143.01135216594;	%7 CCGT
	2	566948.2240065234	226779.28960260932	3	0.0	2892.592979625119	0.0	0.0	58809.397701913185;	%7 OCGT
	2	3883.386907562038	1553.3547630248152	3	0.0	760.1070478688663	0.0	0.0	16189.767503605335;	%7 biomass
	2	0.0	0.0	3	0.0	1.289445522411333	0.0	0.0	21257.35741811583;	%7 onwind
	2	0.0	0.0	3	0.0	1.0174553084742282	0.0	0.0	5803.625312346206;	%7 solar
	2	441815.55484455446	176726.22193782177	3	0.0	2254.160994104869	0.0	0.0	129143.01135216594;	%8 CCGT
	2	579388.2678765142	231755.3071506057	3	0.0	2956.062591206705	0.0	0.0	58809.397701913185;	%8 OCGT
	2	3883.39428963334	1553.357715853336	3	0.0	760.1084927839772	0.0	0.0	16189.767503605335;	%8 biomass
	2	0.0	0.0	3	0.0	1.303441411894657	0.0	0.0	21257.35741811583;	%8 onwind
	2	0.0	0.0	3	0.0	1.016156583095783	0.0	0.0	5803.625312346205;	%8 solar
	2	423663.6586113629	169465.46344454517	3	0.0	2161.5492786294026	0.0	0.0	129143.01135216594;	%9 CCGT
	2	571638.9269332109	228655.57077328436	3	0.0	2916.5251374143413	0.0	0.0	58809.397701913185;	%9 OCGT
	2	3883.2418842745947	1553.296753709838	3	0.0	760.0786620228215	0.0	0.0	16189.767503605335;	%9 biomass
	2	0.0	0.0	3	0.0	1.2676423411133135	0.0	0.0	21257.35741811583;	%9 onwind
	2	0.0	0.0	3	0.0	1.0279996482603462	0.0	0.0	5803.625312346205;	%9 solar
	2	0.0	0.0	3	0.0	0.4930724332539114	0.0	0.0	1032.4892748394757;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.3698043249404336	0.0	0.0	1032.4892748394757;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5371930158609275	0.0	0.0	1032.4892748394757;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.4028947618956956	0.0	0.0	1032.4892748394757;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5605321059279206	0.0	0.0	1032.4892748394757;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.42039907944594046	0.0	0.0	1032.4892748394757;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4681006488720658	0.0	0.0	1032.4892748394757;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.35107548665404936	0.0	0.0	1032.4892748394757;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5569040216547296	0.0	0.0	1032.4892748394757;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.4176780162410472	0.0	0.0	1032.4892748394757;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5357420457804054	0.0	0.0	1032.4892748394757;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.40180653433530406	0.0	0.0	1032.4892748394757;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5099553129168626	0.0	0.0	1032.4892748394757;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.38246648468764693	0.0	0.0	1032.4892748394757;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5176428590896757	0.0	0.0	1032.4892748394757;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.5118166158270085	0.0	0.0	1032.4892748394757;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.3838624618702564	0.0	0.0	1032.4892748394757;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4983942014913984	0.0	0.0	1032.4892748394757;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.3737956511185488	0.0	0.0	1032.4892748394757;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.48923355581325756	0.0	0.0	1032.4892748394757;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.36692516685994314	0.0	0.0	1032.4892748394757;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4685127898125384	0.0	0.0	1032.4892748394757;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.3513845923594038	0.0	0.0	1032.4892748394757;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5167550562934854	0.0	0.0	1032.4892748394757;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.3875662922201141	0.0	0.0	1032.4892748394757;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.47435589210890494	0.0	0.0	1032.4892748394757;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.3557669190816787	0.0	0.0	1032.4892748394757;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.47720971446748195	0.0	0.0	1032.4892748394757;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.3579072858506115	0.0	0.0	1032.4892748394757;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5530253895096996	0.0	0.0	1032.4892748394757;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.4147690421322747	0.0	0.0	1032.4892748394757;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5062518214168596	0.0	0.0	1032.4892748394757;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.37968886606264474	0.0	0.0	1032.4892748394757;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5261513216805729	0.0	0.0	1032.4892748394757;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.39461349126042966	0.0	0.0	1032.4892748394757;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5334111560203068	0.0	0.0	1032.4892748394757;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.40005836701523007	0.0	0.0	1032.4892748394757;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4736060587632618	0.0	0.0	1032.4892748394757;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.551272425174086	0.0	0.0	1032.4892748394757;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.4134543188805645	0.0	0.0	1032.4892748394757;	%DE1 76 PHS (pump mode)
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

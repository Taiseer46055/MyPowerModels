function mpc = mpc_000_020
%MPC_000_020	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 0 	Weight: 15
%	Time step: 20

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 15;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1767.9961840906115	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4644.666062051068	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1514.53692251033	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	1649.0178592401173	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	4603.754040170861	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	5316.593133779842	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	250.92893437740602	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2459.3024650769657	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2499.4100064171034	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	3539.7706692636416	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	-106.2867539024822	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	1481.5873043762913	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4907.819995942936	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3287.5123738153884	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3631.67347042373	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	5414.190675358352	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	1205.2330975140762	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4382.560176633543	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5612.059412857696	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3492.5371848462873	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.480260825116547	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	34.67857238967475	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.441905279831676	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.766392899585048	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	30.880334799442185	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.636577286888535	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	37.63632512420975	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.587076760833934	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	22.87211539447738	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.594512959531869	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	35.39546292661783	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.293109841595998	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	33.55173680817492	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.806048903692147	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	28.48363872328843	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	30.58228024215224	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	33.97264052142161	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	26.017391612975675	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.616038195129358	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	31.32146057959848	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.128794960070355	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.111422896831826	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	27.41513157387534	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	28.008957828979035	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	30.378975120814932	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	28.04858571771533	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	1142.1800119596871	456.87200478387484	3	0.0	223.56234330782678	0.0	0.0	4761.696324589805;	%0 biomass
	2	0.0	0.0	3	0.0	0.3700935349980971	0.0	0.0	5514.795544758967;	%0 offwind
	2	0.0	0.0	3	0.0	0.3614622334272965	0.0	0.0	6252.1639465046555;	%0 onwind
	2	0.0	0.0	3	0.0	0.28588291059927845	0.0	0.0	1706.9486212782956;	%0 solar
	2	133243.5418835763	53297.41675343052	3	0.0	679.8139892019199	0.0	0.0	37983.23863298998;	%1 CCGT
	2	162644.96506822016	65057.98602728807	3	0.0	829.821250348062	0.0	0.0	17296.88167703329;	%1 OCGT
	2	1142.134384147449	456.85375365897966	3	0.0	223.55341243833414	0.0	0.0	4761.696324589804;	%1 biomass
	2	0.0	0.0	3	0.0	0.3687002098784729	0.0	0.0	6252.1639465046555;	%1 onwind
	2	0.0	0.0	3	0.0	0.2990965295940781	0.0	0.0	1706.9486212782956;	%1 solar
	2	116749.15732768828	46699.66293107531	3	0.0	595.6589659575932	0.0	0.0	37983.23863298998;	%10 CCGT
	2	169353.1115782995	67741.2446313198	3	0.0	864.0464876443853	0.0	0.0	17296.88167703329;	%10 OCGT
	2	1142.1553948200915	456.86215792803654	3	0.0	223.557524920746	0.0	0.0	4761.696324589804;	%10 biomass
	2	0.0	0.0	3	0.0	0.37039632812104717	0.0	0.0	5444.367073012044;	%10 offwind
	2	0.0	0.0	3	0.0	0.3719956390429314	0.0	0.0	6252.163946504656;	%10 onwind
	2	0.0	0.0	3	0.0	0.2973341701368303	0.0	0.0	1706.9486212782956;	%10 solar
	2	108938.70188205065	43575.48075282026	3	0.0	555.8097034798502	0.0	0.0	37983.23863298998;	%11 CCGT
	2	161719.57321689036	64687.829286756154	3	0.0	825.0998633514814	0.0	0.0	17296.88167703329;	%11 OCGT
	2	1142.18065465007	456.87226186002795	3	0.0	223.56246910355645	0.0	0.0	4761.696324589804;	%11 biomass
	2	0.0	0.0	3	0.0	0.3797133134764515	0.0	0.0	5413.202087000224;	%11 offwind
	2	0.0	0.0	3	0.0	0.36160215446622107	0.0	0.0	6252.1639465046555;	%11 onwind
	2	0.0	0.0	3	0.0	0.301431038390385	0.0	0.0	1706.9486212782956;	%11 solar
	2	160848.0544727454	64339.22178909816	3	0.0	820.6533391466601	0.0	0.0	17296.88167703329;	%12 OCGT
	2	1142.1589467655972	456.8635787062389	3	0.0	223.5582201537673	0.0	0.0	4761.696324589804;	%12 biomass
	2	0.0	0.0	3	0.0	0.37523929369773984	0.0	0.0	6252.1639465046555;	%12 onwind
	2	0.0	0.0	3	0.0	0.3001516086390752	0.0	0.0	1706.9486212782956;	%12 solar
	2	120927.76313314725	48371.1052532589	3	0.0	616.9783833323838	0.0	0.0	37983.23863298998;	%13 CCGT
	2	168215.83412451347	67286.33364980538	3	0.0	858.244051655681	0.0	0.0	17296.88167703329;	%13 OCGT
	2	1142.1484624425598	456.859384977024	3	0.0	223.55616802555488	0.0	0.0	4761.696324589804;	%13 biomass
	2	0.0	0.0	3	0.0	0.38143858117468993	0.0	0.0	6252.1639465046555;	%13 onwind
	2	0.0	0.0	3	0.0	0.3036546325350222	0.0	0.0	1706.9486212782956;	%13 solar
	2	124615.44002909528	49846.176011638105	3	0.0	635.7930613729351	0.0	0.0	37983.23863298998;	%14 CCGT
	2	168095.4019683798	67238.16078735192	3	0.0	857.6296018794887	0.0	0.0	17296.88167703329;	%14 OCGT
	2	1142.1016474136645	456.84065896546576	3	0.0	223.54700477856025	0.0	0.0	4761.696324589804;	%14 biomass
	2	0.0	0.0	3	0.0	0.3708706553317254	0.0	0.0	5701.061274859127;	%14 offwind
	2	0.0	0.0	3	0.0	0.3672000375325234	0.0	0.0	6252.163946504656;	%14 onwind
	2	0.0	0.0	3	0.0	0.2935479626677962	0.0	0.0	1706.9486212782956;	%14 solar
	2	111915.6525568542	44766.26102274168	3	0.0	570.9982273308888	0.0	0.0	37983.23863298998;	%15 CCGT
	2	156972.09029061403	62788.83611624561	3	0.0	800.8780116868063	0.0	0.0	17296.88167703329;	%15 OCGT
	2	1142.1614827830722	456.86459311322886	3	0.0	223.55871653612687	0.0	0.0	4761.696324589804;	%15 biomass
	2	0.0	0.0	3	0.0	0.3754795413456295	0.0	0.0	5949.62610635095;	%15 offwind
	2	0.0	0.0	3	0.0	0.37168095441927923	0.0	0.0	6252.163946504656;	%15 onwind
	2	0.0	0.0	3	0.0	0.2900913540585632	0.0	0.0	1706.9486212782956;	%15 solar
	2	117150.48861951251	46860.195447805	3	0.0	597.7065745893494	0.0	0.0	37983.23863298998;	%16 CCGT
	2	162092.0449825074	64836.81799300296	3	0.0	827.0002295025887	0.0	0.0	17296.88167703329;	%16 OCGT
	2	1142.1339604346128	456.8535841738452	3	0.0	223.55332950374103	0.0	0.0	4761.696324589805;	%16 biomass
	2	0.0	0.0	3	0.0	0.3710514478905071	0.0	0.0	6252.1639465046555;	%16 onwind
	2	0.0	0.0	3	0.0	0.29838476274183856	0.0	0.0	1706.9486212782958;	%16 solar
	2	120445.97340089161	48178.38936035665	3	0.0	614.5202724535286	0.0	0.0	37983.23863298998;	%17 CCGT
	2	158013.73302029076	63205.4932081163	3	0.0	806.1925154096466	0.0	0.0	17296.88167703329;	%17 OCGT
	2	1142.146858266563	456.8587433066252	3	0.0	223.55585403534212	0.0	0.0	4761.696324589804;	%17 biomass
	2	0.0	0.0	3	0.0	0.38238170961621093	0.0	0.0	6252.163946504656;	%17 onwind
	2	0.0	0.0	3	0.0	0.30575987198392285	0.0	0.0	1706.9486212782958;	%17 solar
	2	126254.52569336987	50501.810277347955	3	0.0	644.1557433335197	0.0	0.0	37983.23863298998;	%18 CCGT
	2	159243.4918009606	63697.39672038424	3	0.0	812.4667949028601	0.0	0.0	17296.88167703329;	%18 OCGT
	2	1142.1675310175042	456.8670124070017	3	0.0	223.55990037531888	0.0	0.0	4761.696324589804;	%18 biomass
	2	0.0	0.0	3	0.0	0.3736439292893007	0.0	0.0	6252.1639465046555;	%18 onwind
	2	0.0	0.0	3	0.0	0.2999915134545782	0.0	0.0	1706.9486212782956;	%18 solar
	2	118302.2517906427	47320.90071625708	3	0.0	603.5829172991974	0.0	0.0	37983.23863298998;	%19 CCGT
	2	165752.14695737176	66300.85878294871	3	0.0	845.674219170264	0.0	0.0	17296.88167703329;	%19 OCGT
	2	1142.1412132550709	456.8564853020284	3	0.0	223.55474912019395	0.0	0.0	4761.696324589804;	%19 biomass
	2	0.0	0.0	3	0.0	0.3812105379220506	0.0	0.0	6252.1639465046555;	%19 onwind
	2	0.0	0.0	3	0.0	0.2981874573802173	0.0	0.0	1706.9486212782956;	%19 solar
	2	126653.72897775119	50661.49159110049	3	0.0	646.1924947844449	0.0	0.0	37983.23863298998;	%2 CCGT
	2	167388.69711868904	66955.47884747562	3	0.0	854.0239648912706	0.0	0.0	17296.88167703329;	%2 OCGT
	2	1142.1356275136502	456.85425100546007	3	0.0	223.55365580615583	0.0	0.0	4761.696324589804;	%2 biomass
	2	0.0	0.0	3	0.0	0.37568046672533634	0.0	0.0	6252.1639465046555;	%2 onwind
	2	0.0	0.0	3	0.0	0.299145274773518	0.0	0.0	1706.9486212782956;	%2 solar
	2	138118.82333094106	55247.52933237643	3	0.0	704.6878741374544	0.0	0.0	37983.23863298998;	%3 CCGT
	2	171809.51708645458	68723.80683458183	3	0.0	876.5791688084416	0.0	0.0	17296.88167703329;	%3 OCGT
	2	1142.1022637072722	456.8409054829089	3	0.0	223.54712540756944	0.0	0.0	4761.696324589804;	%3 biomass
	2	0.0	0.0	3	0.0	0.3716146064736323	0.0	0.0	5923.136552979391;	%3 offwind
	2	0.0	0.0	3	0.0	0.3823037559488863	0.0	0.0	6252.163946504656;	%3 onwind
	2	0.0	0.0	3	0.0	0.299552285917403	0.0	0.0	1706.9486212782956;	%3 solar
	2	120399.97954928495	48159.991819713985	3	0.0	614.2856099453314	0.0	0.0	37983.23863298998;	%4 CCGT
	2	164657.24777804303	65862.89911121721	3	0.0	840.0879988675665	0.0	0.0	17296.88167703329;	%4 OCGT
	2	1142.1701001926344	456.8680400770538	3	0.0	223.56040324772647	0.0	0.0	4761.696324589804;	%4 biomass
	2	0.0	0.0	3	0.0	0.3739471144007091	0.0	0.0	6252.1639465046555;	%4 onwind
	2	0.0	0.0	3	0.0	0.2931680798200606	0.0	0.0	1706.9486212782956;	%4 solar
	2	123666.47064391641	49466.58825756657	3	0.0	630.9513808363082	0.0	0.0	37983.23863298998;	%5 CCGT
	2	161690.29879418484	64676.11951767393	3	0.0	824.9505040519633	0.0	0.0	17296.88167703329;	%5 OCGT
	2	1142.154903753824	456.86196150152966	3	0.0	223.55742880286243	0.0	0.0	4761.696324589804;	%5 biomass
	2	0.0	0.0	3	0.0	0.36940091066952363	0.0	0.0	6252.1639465046555;	%5 onwind
	2	0.0	0.0	3	0.0	0.3002550818687922	0.0	0.0	1706.9486212782956;	%5 solar
	2	121253.76061771157	48501.50424708463	3	0.0	618.6416358046508	0.0	0.0	37983.23863298998;	%6 CCGT
	2	167465.42729470338	66986.17091788136	3	0.0	854.4154453811398	0.0	0.0	17296.88167703329;	%6 OCGT
	2	1142.1738429847942	456.86953719391767	3	0.0	223.5611358357397	0.0	0.0	4761.696324589804;	%6 biomass
	2	0.0	0.0	3	0.0	0.3751453695869937	0.0	0.0	6252.1639465046555;	%6 onwind
	2	0.0	0.0	3	0.0	0.3062828535592847	0.0	0.0	1706.9486212782956;	%6 solar
	2	114931.01077314315	45972.40430925726	3	0.0	586.3827080262406	0.0	0.0	37983.23863298998;	%7 CCGT
	2	166749.47764897745	66699.79105959098	3	0.0	850.7626410662115	0.0	0.0	17296.88167703329;	%7 OCGT
	2	1142.1726198711876	456.86904794847504	3	0.0	223.5608964320195	0.0	0.0	4761.696324589804;	%7 biomass
	2	0.0	0.0	3	0.0	0.37924868306215676	0.0	0.0	6252.1639465046555;	%7 onwind
	2	0.0	0.0	3	0.0	0.2992515613159495	0.0	0.0	1706.9486212782958;	%7 solar
	2	129945.75142486895	51978.300569947576	3	0.0	662.9885276779028	0.0	0.0	37983.23863298998;	%8 CCGT
	2	170408.3140813277	68163.32563253108	3	0.0	869.430173884325	0.0	0.0	17296.88167703329;	%8 OCGT
	2	1142.1747910686295	456.86991642745176	3	0.0	223.56132140705213	0.0	0.0	4761.696324589804;	%8 biomass
	2	0.0	0.0	3	0.0	0.3833651211454873	0.0	0.0	6252.1639465046555;	%8 onwind
	2	0.0	0.0	3	0.0	0.29886958326346563	0.0	0.0	1706.9486212782956;	%8 solar
	2	124606.95841510674	49842.7833660427	3	0.0	635.7497878321772	0.0	0.0	37983.23863298998;	%9 CCGT
	2	168129.09615682674	67251.6384627307	3	0.0	857.801511004218	0.0	0.0	17296.88167703329;	%9 OCGT
	2	1142.129965963116	456.8519863852465	3	0.0	223.55254765377103	0.0	0.0	4761.696324589804;	%9 biomass
	2	0.0	0.0	3	0.0	0.37283598268038637	0.0	0.0	6252.1639465046555;	%9 onwind
	2	0.0	0.0	3	0.0	0.30235283772363125	0.0	0.0	1706.9486212782956;	%9 solar
	2	0.0	0.0	3	0.0	0.14502130389820925	0.0	0.0	303.67331612925756;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.10876597792365694	0.0	0.0	303.67331612925756;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15799794584144924	0.0	0.0	303.67331612925756;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.11849845938108694	0.0	0.0	303.67331612925756;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16486238409644724	0.0	0.0	303.67331612925756;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.12364678807233542	0.0	0.0	303.67331612925756;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13767666143296053	0.0	0.0	303.67331612925756;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.1032574960747204	0.0	0.0	303.67331612925756;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1637953004866852	0.0	0.0	303.67331612925756;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.12284647536501389	0.0	0.0	303.67331612925756;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15757118993541333	0.0	0.0	303.67331612925756;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.11817839245156	0.0	0.0	303.67331612925756;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1499868567402537	0.0	0.0	303.67331612925756;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.11249014255519027	0.0	0.0	303.67331612925756;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15224789973225755	0.0	0.0	303.67331612925756;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.15053429877264957	0.0	0.0	303.67331612925756;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.11290072407948717	0.0	0.0	303.67331612925756;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1465865298504113	0.0	0.0	303.67331612925756;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.10993989738780846	0.0	0.0	303.67331612925756;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.14389222229801693	0.0	0.0	303.67331612925756;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.1079191667235127	0.0	0.0	303.67331612925756;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13779787935662896	0.0	0.0	303.67331612925756;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.10334840951747172	0.0	0.0	303.67331612925756;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15198678126278983	0.0	0.0	303.67331612925756;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.11399008594709237	0.0	0.0	303.67331612925756;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13951643885556028	0.0	0.0	303.67331612925756;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.1046373291416702	0.0	0.0	303.67331612925756;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1403557983727888	0.0	0.0	303.67331612925756;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.10526684877959161	0.0	0.0	303.67331612925756;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.16265452632638225	0.0	0.0	303.67331612925756;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.1219908947447867	0.0	0.0	303.67331612925756;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1488975945343705	0.0	0.0	303.67331612925756;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.11167319590077787	0.0	0.0	303.67331612925756;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15475038872958025	0.0	0.0	303.67331612925756;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.11606279154718518	0.0	0.0	303.67331612925756;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.15688563412361964	0.0	0.0	303.67331612925756;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.11766422559271472	0.0	0.0	303.67331612925756;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13929589963625347	0.0	0.0	303.67331612925756;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.16213894858061353	0.0	0.0	303.67331612925756;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.12160421143546016	0.0	0.0	303.67331612925756;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_002_003
%MPC_002_003	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 2 	Weight: 29
%	Time step: 3

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 29;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	191.52571718117667	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	3385.625105743039	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	512.5452036237334	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-362.2885961341632	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	1224.125403381565	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	4316.448387850722	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	117.9295487857513	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	1765.4436586033396	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	1195.6756236557987	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	896.150895039687	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	2309.9752129483686	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	2095.5671500970534	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	3819.5777552543473	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	1582.2968081070762	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	1247.7991020102738	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-2728.524439545151	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	1438.20209528718	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	2988.7420258392535	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3020.2276822396702	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	2567.0568176392294	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.2541282108043057	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.9278890686715657	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.5531739626961139	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7717078511906825	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.243200220464862	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.570591194858698	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.555511921112217	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.204636730472167	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.7267163023330784	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.9071797954728142	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.4421396993198154	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.925631369269709	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5848341385850407	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.317244307527674	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.834200549966865	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.508738079349389	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6664391447419757	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.655853626101631	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.885959270705571	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.307453476184432	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.788511407493525	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.438426054192182	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.2978752173841452	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.8409817516403786	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.494077058929699	0.0	0.0	133	133	656	38.0	2;	%9 onwind
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
	2	2208.2146897887283	883.2858759154914	3	0.0	432.2205303951318	0.0	0.0	9205.946227540291;	%0 biomass
	2	0.0	0.0	3	0.0	0.7155141676629877	0.0	0.0	10661.938053200669;	%0 offwind
	2	0.0	0.0	3	0.0	0.6988269846261066	0.0	0.0	12087.516963242335;	%0 onwind
	2	0.0	0.0	3	0.0	0.5527069604919383	0.0	0.0	3300.1006678047047;	%0 solar
	2	257604.18097491423	103041.67238996568	3	0.0	1314.3070457903784	0.0	0.0	73434.26135711397;	%1 CCGT
	2	314446.9324652257	125778.77298609026	3	0.0	1604.3210840062532	0.0	0.0	33440.63790893102;	%1 OCGT
	2	2208.1264760184017	883.2505904073606	3	0.0	432.203264047446	0.0	0.0	9205.94622754029;	%1 biomass
	2	0.0	0.0	3	0.0	0.7128204057650476	0.0	0.0	12087.516963242335;	%1 onwind
	2	0.0	0.0	3	0.0	0.578253290548551	0.0	0.0	3300.1006678047047;	%1 solar
	2	225715.03750019736	90286.01500007894	3	0.0	1151.6073341846802	0.0	0.0	73434.26135711397;	%10 CCGT
	2	327416.0157180457	130966.40628721828	3	0.0	1670.489876112478	0.0	0.0	33440.63790893102;	%10 OCGT
	2	2208.167096652177	883.2668386608707	3	0.0	432.2112148467756	0.0	0.0	9205.94622754029;	%10 biomass
	2	0.0	0.0	3	0.0	0.7160995677006912	0.0	0.0	10525.776341156617;	%10 offwind
	2	0.0	0.0	3	0.0	0.7191915688163342	0.0	0.0	12087.516963242335;	%10 onwind
	2	0.0	0.0	3	0.0	0.5748460622645385	0.0	0.0	3300.1006678047047;	%10 solar
	2	210614.82363863126	84245.9294554525	3	0.0	1074.5654267277105	0.0	0.0	73434.26135711397;	%11 CCGT
	2	312657.84155265475	125063.13662106189	3	0.0	1595.1930691461973	0.0	0.0	33440.63790893102;	%11 OCGT
	2	2208.2159323234687	883.2863729293873	3	0.0	432.2207736002091	0.0	0.0	9205.94622754029;	%11 biomass
	2	0.0	0.0	3	0.0	0.7341124060544729	0.0	0.0	10465.5240348671;	%11 offwind
	2	0.0	0.0	3	0.0	0.6990974986346941	0.0	0.0	12087.516963242335;	%11 onwind
	2	0.0	0.0	3	0.0	0.5827666742214109	0.0	0.0	3300.1006678047047;	%11 solar
	2	310972.9053139745	124389.16212558978	3	0.0	1586.596455683543	0.0	0.0	33440.63790893102;	%12 OCGT
	2	2208.1739637468213	883.2695854987285	3	0.0	432.21255896395013	0.0	0.0	9205.94622754029;	%12 biomass
	2	0.0	0.0	3	0.0	0.725462634482297	0.0	0.0	12087.516963242335;	%12 onwind
	2	0.0	0.0	3	0.0	0.5802931100355453	0.0	0.0	3300.1006678047047;	%12 solar
	2	233793.67539075136	93517.47015630055	3	0.0	1192.824874442609	0.0	0.0	73434.26135711397;	%13 CCGT
	2	325217.2793073927	130086.91172295708	3	0.0	1659.271833200983	0.0	0.0	33440.63790893102;	%13 OCGT
	2	2208.1536940556157	883.2614776222463	3	0.0	432.2085915160728	0.0	0.0	9205.94622754029;	%13 biomass
	2	0.0	0.0	3	0.0	0.7374479236044005	0.0	0.0	12087.516963242335;	%13 onwind
	2	0.0	0.0	3	0.0	0.5870656229010429	0.0	0.0	3300.1006678047047;	%13 solar
	2	240923.18405625084	96369.27362250035	3	0.0	1229.199918654341	0.0	0.0	73434.26135711397;	%14 CCGT
	2	324984.4438055343	129993.7775222137	3	0.0	1658.0838969670115	0.0	0.0	33440.63790893102;	%14 OCGT
	2	2208.0631849997512	883.2252739999004	3	0.0	432.1908759052165	0.0	0.0	9205.94622754029;	%14 biomass
	2	0.0	0.0	3	0.0	0.7170166003080025	0.0	0.0	11022.051798060977;	%14 offwind
	2	0.0	0.0	3	0.0	0.7099200725628785	0.0	0.0	12087.516963242335;	%14 onwind
	2	0.0	0.0	3	0.0	0.5675260611577393	0.0	0.0	3300.1006678047047;	%14 solar
	2	216370.26160991812	86548.10464396725	3	0.0	1103.9299061730514	0.0	0.0	73434.26135711397;	%15 CCGT
	2	303479.37456185377	121391.74982474152	3	0.0	1548.3641559278255	0.0	0.0	33440.63790893102;	%15 OCGT
	2	2208.1788667139394	883.2715466855758	3	0.0	432.21351863651194	0.0	0.0	9205.94622754029;	%15 biomass
	2	0.0	0.0	3	0.0	0.725927113268217	0.0	0.0	11502.610472278504;	%15 offwind
	2	0.0	0.0	3	0.0	0.7185831785439398	0.0	0.0	12087.516963242335;	%15 onwind
	2	0.0	0.0	3	0.0	0.5608432845132222	0.0	0.0	3300.1006678047047;	%15 solar
	2	226490.94466439087	90596.37786575634	3	0.0	1155.5660442060757	0.0	0.0	73434.26135711397;	%16 CCGT
	2	313377.95363284764	125351.18145313904	3	0.0	1598.8671103716715	0.0	0.0	33440.63790893102;	%16 OCGT
	2	2208.125656840252	883.2502627361007	3	0.0	432.20310370723263	0.0	0.0	9205.946227540291;	%16 biomass
	2	0.0	0.0	3	0.0	0.7173661325883137	0.0	0.0	12087.516963242335;	%16 onwind
	2	0.0	0.0	3	0.0	0.5768772079675546	0.0	0.0	3300.100667804705;	%16 solar
	2	232862.21524172378	93144.88609668952	3	0.0	1188.0725267434887	0.0	0.0	73434.26135711397;	%17 CCGT
	2	305493.21717256214	122197.28686902484	3	0.0	1558.6388631253167	0.0	0.0	33440.63790893102;	%17 OCGT
	2	2208.1505926486884	883.2602370594753	3	0.0	432.2079844683281	0.0	0.0	9205.94622754029;	%17 biomass
	2	0.0	0.0	3	0.0	0.7392713052580079	0.0	0.0	12087.516963242335;	%17 onwind
	2	0.0	0.0	3	0.0	0.5911357525022508	0.0	0.0	3300.100667804705;	%17 solar
	2	244092.08300718176	97636.8332028727	3	0.0	1245.3677704448048	0.0	0.0	73434.26135711397;	%18 CCGT
	2	307870.7508151905	123148.30032607619	3	0.0	1570.7691368121962	0.0	0.0	33440.63790893102;	%18 OCGT
	2	2208.190559967175	883.2762239868699	3	0.0	432.2158073922832	0.0	0.0	9205.94622754029;	%18 biomass
	2	0.0	0.0	3	0.0	0.722378263292648	0.0	0.0	12087.516963242335;	%18 onwind
	2	0.0	0.0	3	0.0	0.5799835926788511	0.0	0.0	3300.1006678047047;	%18 solar
	2	228717.68679524254	91487.07471809702	3	0.0	1166.926973445115	0.0	0.0	73434.26135711397;	%19 CCGT
	2	320454.1507842521	128181.66031370084	3	0.0	1634.9701570625105	0.0	0.0	33440.63790893102;	%19 OCGT
	2	2208.139678959804	883.2558715839215	3	0.0	432.20584829904163	0.0	0.0	9205.94622754029;	%19 biomass
	2	0.0	0.0	3	0.0	0.7370070399826313	0.0	0.0	12087.516963242335;	%19 onwind
	2	0.0	0.0	3	0.0	0.5764957509350868	0.0	0.0	3300.1006678047047;	%19 solar
	2	244863.87602365232	97945.55040946094	3	0.0	1249.3054899165934	0.0	0.0	73434.26135711397;	%2 CCGT
	2	323618.1477627988	129447.25910511953	3	0.0	1651.1129987897898	0.0	0.0	33440.63790893102;	%2 OCGT
	2	2208.1288798597234	883.2515519438895	3	0.0	432.20373455856793	0.0	0.0	9205.94622754029;	%2 biomass
	2	0.0	0.0	3	0.0	0.7263155690023169	0.0	0.0	12087.516963242335;	%2 onwind
	2	0.0	0.0	3	0.0	0.5783475312288014	0.0	0.0	3300.1006678047047;	%2 solar
	2	267029.72510648606	106811.89004259443	3	0.0	1362.3965566657453	0.0	0.0	73434.26135711397;	%3 CCGT
	2	332165.0663671455	132866.02654685822	3	0.0	1694.7197263629873	0.0	0.0	33440.63790893102;	%3 OCGT
	2	2208.064376500726	883.2257506002906	3	0.0	432.19110912130094	0.0	0.0	9205.94622754029;	%3 biomass
	2	0.0	0.0	3	0.0	0.7184549058490224	0.0	0.0	11451.397335760155;	%3 offwind
	2	0.0	0.0	3	0.0	0.7391205948345135	0.0	0.0	12087.516963242335;	%3 onwind
	2	0.0	0.0	3	0.0	0.5791344194403124	0.0	0.0	3300.1006678047047;	%3 solar
	2	232773.29379528426	93109.3175181137	3	0.0	1187.6188458943075	0.0	0.0	73434.26135711397;	%4 CCGT
	2	318337.3457042165	127334.93828168661	3	0.0	1624.170131143962	0.0	0.0	33440.63790893102;	%4 OCGT
	2	2208.1955270390936	883.2782108156373	3	0.0	432.2167796122712	0.0	0.0	9205.94622754029;	%4 biomass
	2	0.0	0.0	3	0.0	0.7229644211747043	0.0	0.0	12087.516963242335;	%4 onwind
	2	0.0	0.0	3	0.0	0.5667916209854504	0.0	0.0	3300.1006678047047;	%4 solar
	2	239088.50991157175	95635.4039646287	3	0.0	1219.8393362835293	0.0	0.0	73434.26135711397;	%5 CCGT
	2	312601.244335424	125040.49773416959	3	0.0	1594.9043078337957	0.0	0.0	33440.63790893102;	%5 OCGT
	2	2208.1661472573933	883.2664589029573	3	0.0	432.21102901886735	0.0	0.0	9205.94622754029;	%5 biomass
	2	0.0	0.0	3	0.0	0.714175093961079	0.0	0.0	12087.516963242335;	%5 onwind
	2	0.0	0.0	3	0.0	0.5804931582796649	0.0	0.0	3300.1006678047047;	%5 solar
	2	234423.93719424237	93769.57487769696	3	0.0	1196.0404958889917	0.0	0.0	73434.26135711397;	%6 CCGT
	2	323766.4927697599	129506.59710790394	3	0.0	1651.8698610702036	0.0	0.0	33440.63790893102;	%6 OCGT
	2	2208.2027631039355	883.2811052415741	3	0.0	432.21819594909675	0.0	0.0	9205.94622754029;	%6 biomass
	2	0.0	0.0	3	0.0	0.7252810478681877	0.0	0.0	12087.516963242335;	%6 onwind
	2	0.0	0.0	3	0.0	0.592146850214617	0.0	0.0	3300.1006678047047;	%6 solar
	2	222199.95416141007	88879.98166456404	3	0.0	1133.6732355173983	0.0	0.0	73434.26135711397;	%7 CCGT
	2	322382.32345468976	128952.92938187589	3	0.0	1644.8077727280088	0.0	0.0	33440.63790893102;	%7 OCGT
	2	2208.2003984176295	883.2801593670517	3	0.0	432.21773310190434	0.0	0.0	9205.94622754029;	%7 biomass
	2	0.0	0.0	3	0.0	0.7332141205868364	0.0	0.0	12087.516963242335;	%7 onwind
	2	0.0	0.0	3	0.0	0.578553018544169	0.0	0.0	3300.100667804705;	%7 solar
	2	251228.45275474666	100491.38110189865	3	0.0	1281.7778201772785	0.0	0.0	73434.26135711397;	%8 CCGT
	2	329456.0738905669	131782.42955622677	3	0.0	1680.8983361763617	0.0	0.0	33440.63790893102;	%8 OCGT
	2	2208.2045960660166	883.2818384264067	3	0.0	432.2185547203008	0.0	0.0	9205.94622754029;	%8 biomass
	2	0.0	0.0	3	0.0	0.7411725675479421	0.0	0.0	12087.516963242335;	%8 onwind
	2	0.0	0.0	3	0.0	0.5778145276427001	0.0	0.0	3300.1006678047047;	%8 solar
	2	240906.78626920635	96362.71450768254	3	0.0	1229.1162564755425	0.0	0.0	73434.26135711397;	%9 CCGT
	2	325049.58590319834	130019.83436127934	3	0.0	1658.4162546081548	0.0	0.0	33440.63790893102;	%9 OCGT
	2	2208.117934195358	883.2471736781432	3	0.0	432.20159213062396	0.0	0.0	9205.94622754029;	%9 biomass
	2	0.0	0.0	3	0.0	0.7208162331820803	0.0	0.0	12087.516963242335;	%9 onwind
	2	0.0	0.0	3	0.0	0.5845488195990204	0.0	0.0	3300.1006678047047;	%9 solar
	2	0.0	0.0	3	0.0	0.2803745208698712	0.0	0.0	587.1017445165646;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.2102808906524034	0.0	0.0	587.1017445165646;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30546269529346853	0.0	0.0	587.1017445165646;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.22909702147010141	0.0	0.0	587.1017445165646;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31873394258646465	0.0	0.0	587.1017445165646;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.2390504569398485	0.0	0.0	587.1017445165646;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26617487877039037	0.0	0.0	587.1017445165646;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.19963115907779277	0.0	0.0	587.1017445165646;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31667091427425803	0.0	0.0	587.1017445165646;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.2375031857056935	0.0	0.0	587.1017445165646;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30463763387513243	0.0	0.0	587.1017445165646;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.22847822540634932	0.0	0.0	587.1017445165646;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2899745896978238	0.0	0.0	587.1017445165646;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.21748094227336784	0.0	0.0	587.1017445165646;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2943459394823646	0.0	0.0	587.1017445165646;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.2910329776271225	0.0	0.0	587.1017445165646;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.21827473322034185	0.0	0.0	587.1017445165646;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28340062437746183	0.0	0.0	587.1017445165646;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.21255046828309637	0.0	0.0	587.1017445165646;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.27819162977616607	0.0	0.0	587.1017445165646;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.20864372233212455	0.0	0.0	587.1017445165646;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26640923342281597	0.0	0.0	587.1017445165646;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.19980692506711198	0.0	0.0	587.1017445165646;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29384111044139366	0.0	0.0	587.1017445165646;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.22038083283104526	0.0	0.0	587.1017445165646;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2697317817874165	0.0	0.0	587.1017445165646;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.2022988363405624	0.0	0.0	587.1017445165646;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.271354543520725	0.0	0.0	587.1017445165646;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.20351590764054378	0.0	0.0	587.1017445165646;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.314465417564339	0.0	0.0	587.1017445165646;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.23584906317325424	0.0	0.0	587.1017445165646;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2878686827664496	0.0	0.0	587.1017445165646;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.2159015120748372	0.0	0.0	587.1017445165646;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29918408487718845	0.0	0.0	587.1017445165646;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.22438806365789132	0.0	0.0	587.1017445165646;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3033122259723313	0.0	0.0	587.1017445165646;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.2274841694792485	0.0	0.0	587.1017445165646;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2693054059634234	0.0	0.0	587.1017445165646;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.3134686339225195	0.0	0.0	587.1017445165646;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.23510147544188964	0.0	0.0	587.1017445165646;	%DE1 76 PHS (pump mode)
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

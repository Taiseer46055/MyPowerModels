function mpc = mpc_008_009
%MPC_008_009	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 8 	Weight: 50
%	Time step: 9

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 50;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	400.99707343860456	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5675.581630973153	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	2128.9341592881174	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	1251.360050644676	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	-1129.123276649642	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	7435.315157216027	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	391.2672808754173	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3206.83708110525	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	3056.107389729246	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	5136.830399921481	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	37.949323808283424	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	6378.889057106292	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6394.108821047033	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4408.418106806929	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2685.7722309947044	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	1868.2352590402447	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8864.067828597797	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	6057.7320196207975	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4090.249253274075	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4293.762656467583	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8254988539685858	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.3166687768207892	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.169332142041455	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.901772389954276	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5132393516065865	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.785398681964841	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.705728574590541	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.1944933286601739	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.0899543454241711	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.09386530140955	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.807052772082478	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.4883319077079	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.4948212747422512	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6964575310810076	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.519018023958298	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8618450457474747	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.9306027427493744	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.925365115030186	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.970926106649214	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.058789119773426525	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.688917122455647	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.7707214243914526	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.561615735027742	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.760099308387643	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.919296497733709	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.105831930452464	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.184245410936353	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5805604885263193	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.02357665201547	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.94112298869461	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.33853743248193224	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.580077641077782	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.14836137534711	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.640894150551351	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.457761473021492	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.476271409413802	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.209388704744796	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.30297874849727074	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.632977153201718	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	3807.2667065322903	1522.906682612916	3	0.0	745.2078110260893	0.0	0.0	15872.321081966018;	%0 biomass
	2	0.0	0.0	3	0.0	1.2336451166603237	0.0	0.0	18382.65181586322;	%0 offwind
	2	0.0	0.0	3	0.0	1.2048741114243218	0.0	0.0	20840.546488348853;	%0 onwind
	2	0.0	0.0	3	0.0	0.9529430353309282	0.0	0.0	5689.828737594318;	%0 solar
	2	444145.13961192104	177658.05584476842	3	0.0	2266.0466306730664	0.0	0.0	126610.79544329994;	%1 CCGT
	2	542149.8835607339	216859.95342429355	3	0.0	2766.07083449354	0.0	0.0	57656.27225677763;	%1 OCGT
	2	3807.11461382483	1522.8458455299321	3	0.0	745.1780414611138	0.0	0.0	15872.321081966016;	%1 biomass
	2	0.0	0.0	3	0.0	1.2290006995949097	0.0	0.0	20840.546488348853;	%1 onwind
	2	0.0	0.0	3	0.0	0.9969884319802604	0.0	0.0	5689.828737594318;	%1 solar
	2	389163.857758961	155665.54310358438	3	0.0	1985.5298865253108	0.0	0.0	126610.79544329994;	%10 CCGT
	2	564510.371927665	225804.14877106601	3	0.0	2880.154958814617	0.0	0.0	57656.27225677763;	%10 OCGT
	2	3807.184649400305	1522.8738597601218	3	0.0	745.1917497358199	0.0	0.0	15872.321081966016;	%10 biomass
	2	0.0	0.0	3	0.0	1.2346544270701572	0.0	0.0	18147.89024337348;	%10 offwind
	2	0.0	0.0	3	0.0	1.2399854634764382	0.0	0.0	20840.546488348853;	%10 onwind
	2	0.0	0.0	3	0.0	0.9911139004561009	0.0	0.0	5689.828737594318;	%10 solar
	2	363129.0062735022	145251.60250940087	3	0.0	1852.6990115995009	0.0	0.0	126610.79544329994;	%11 CCGT
	2	539065.2440563013	215626.0976225205	3	0.0	2750.332877838271	0.0	0.0	57656.27225677763;	%11 OCGT
	2	3807.2688488335666	1522.9075395334264	3	0.0	745.2082303451881	0.0	0.0	15872.321081966016;	%11 biomass
	2	0.0	0.0	3	0.0	1.265711044921505	0.0	0.0	18044.006956667414;	%11 offwind
	2	0.0	0.0	3	0.0	1.2053405148874037	0.0	0.0	20840.546488348853;	%11 onwind
	2	0.0	0.0	3	0.0	1.00477012796795	0.0	0.0	5689.828737594318;	%11 solar
	2	536160.181575818	214464.0726303272	3	0.0	2735.511130488867	0.0	0.0	57656.27225677763;	%12 OCGT
	2	3807.196489218657	1522.8785956874628	3	0.0	745.1940671792244	0.0	0.0	15872.321081966016;	%12 biomass
	2	0.0	0.0	3	0.0	1.2507976456591328	0.0	0.0	20840.546488348853;	%12 onwind
	2	0.0	0.0	3	0.0	1.0005053621302507	0.0	0.0	5689.828737594318;	%12 solar
	2	403092.5437771575	161237.01751086302	3	0.0	2056.594611107946	0.0	0.0	126610.79544329994;	%13 CCGT
	2	560719.4470817116	224287.77883268462	3	0.0	2860.8135055189364	0.0	0.0	57656.27225677763;	%13 OCGT
	2	3807.1615414751996	1522.8646165900798	3	0.0	745.1872267518496	0.0	0.0	15872.321081966016;	%13 biomass
	2	0.0	0.0	3	0.0	1.2714619372489664	0.0	0.0	20840.546488348853;	%13 onwind
	2	0.0	0.0	3	0.0	1.012182108450074	0.0	0.0	5689.828737594318;	%13 solar
	2	415384.80009698425	166153.9200387937	3	0.0	2119.31020457645	0.0	0.0	126610.79544329994;	%14 CCGT
	2	560318.006561266	224127.20262450638	3	0.0	2858.7653395982957	0.0	0.0	57656.27225677763;	%14 OCGT
	2	3807.0054913788817	1522.8021965515525	3	0.0	745.1566825952009	0.0	0.0	15872.321081966016;	%14 biomass
	2	0.0	0.0	3	0.0	1.236235517772418	0.0	0.0	19003.537582863755;	%14 offwind
	2	0.0	0.0	3	0.0	1.2240001251084112	0.0	0.0	20840.546488348853;	%14 onwind
	2	0.0	0.0	3	0.0	0.9784932088926539	0.0	0.0	5689.828737594318;	%14 solar
	2	373052.175189514	149220.8700758056	3	0.0	1903.3274244362956	0.0	0.0	126610.79544329994;	%15 CCGT
	2	523240.30096871336	209296.12038748537	3	0.0	2669.593372289354	0.0	0.0	57656.27225677763;	%15 OCGT
	2	3807.2049426102403	1522.8819770440962	3	0.0	745.1957217870896	0.0	0.0	15872.321081966016;	%15 biomass
	2	0.0	0.0	3	0.0	1.2515984711520984	0.0	0.0	19832.087021169835;	%15 offwind
	2	0.0	0.0	3	0.0	1.2389365147309308	0.0	0.0	20840.546488348853;	%15 onwind
	2	0.0	0.0	3	0.0	0.9669711801952107	0.0	0.0	5689.828737594318;	%15 solar
	2	390501.6287317084	156200.65149268333	3	0.0	1992.355248631165	0.0	0.0	126610.79544329994;	%16 CCGT
	2	540306.8166083579	216122.7266433432	3	0.0	2756.6674316752956	0.0	0.0	57656.27225677763;	%16 OCGT
	2	3807.11320144871	1522.8452805794839	3	0.0	745.1777650124701	0.0	0.0	15872.321081966018;	%16 biomass
	2	0.0	0.0	3	0.0	1.2368381596350235	0.0	0.0	20840.546488348853;	%16 onwind
	2	0.0	0.0	3	0.0	0.9946158758061285	0.0	0.0	5689.828737594319;	%16 solar
	2	401486.57800297206	160594.63120118884	3	0.0	2048.400908178429	0.0	0.0	126610.79544329994;	%17 CCGT
	2	526712.4434009691	210684.97736038765	3	0.0	2687.308384698822	0.0	0.0	57656.27225677763;	%17 OCGT
	2	3807.1561942218764	1522.8624776887505	3	0.0	745.1861801178071	0.0	0.0	15872.321081966016;	%17 biomass
	2	0.0	0.0	3	0.0	1.2746056987207033	0.0	0.0	20840.546488348853;	%17 onwind
	2	0.0	0.0	3	0.0	1.0191995732797428	0.0	0.0	5689.828737594319;	%17 solar
	2	420848.4189778996	168339.36759115983	3	0.0	2147.1858111117326	0.0	0.0	126610.79544329994;	%18 CCGT
	2	530811.6393365354	212324.65573461412	3	0.0	2708.2226496762005	0.0	0.0	57656.27225677763;	%18 OCGT
	2	3807.2251033916805	1522.8900413566723	3	0.0	745.1996679177296	0.0	0.0	15872.321081966016;	%18 biomass
	2	0.0	0.0	3	0.0	1.245479764297669	0.0	0.0	20840.546488348853;	%18 onwind
	2	0.0	0.0	3	0.0	0.9999717115152607	0.0	0.0	5689.828737594318;	%18 solar
	2	394340.8393021423	157736.33572085694	3	0.0	2011.9430576639913	0.0	0.0	126610.79544329994;	%19 CCGT
	2	552507.1565245725	221002.86260982903	3	0.0	2818.91406390088	0.0	0.0	57656.27225677763;	%19 OCGT
	2	3807.1373775169027	1522.8549510067612	3	0.0	745.1824970673132	0.0	0.0	15872.321081966016;	%19 biomass
	2	0.0	0.0	3	0.0	1.270701793073502	0.0	0.0	20840.546488348853;	%19 onwind
	2	0.0	0.0	3	0.0	0.9939581912673909	0.0	0.0	5689.828737594318;	%19 solar
	2	422179.096592504	168871.63863700163	3	0.0	2153.9749826148163	0.0	0.0	126610.79544329994;	%2 CCGT
	2	557962.3237289635	223184.9294915854	3	0.0	2846.7465496375685	0.0	0.0	57656.27225677763;	%2 OCGT
	2	3807.1187583788337	1522.8475033515335	3	0.0	745.1788526871861	0.0	0.0	15872.321081966016;	%2 biomass
	2	0.0	0.0	3	0.0	1.2522682224177877	0.0	0.0	20840.546488348853;	%2 onwind
	2	0.0	0.0	3	0.0	0.9971509159117267	0.0	0.0	5689.828737594318;	%2 solar
	2	460396.0777698036	184158.43110792144	3	0.0	2348.9595804581813	0.0	0.0	126610.79544329994;	%3 CCGT
	2	572698.390288182	229079.35611527276	3	0.0	2921.9305626948058	0.0	0.0	57656.27225677763;	%3 OCGT
	2	3807.007545690907	1522.803018276363	3	0.0	745.1570846918981	0.0	0.0	15872.321081966016;	%3 biomass
	2	0.0	0.0	3	0.0	1.2387153549121077	0.0	0.0	19743.788509931303;	%3 offwind
	2	0.0	0.0	3	0.0	1.2743458531629543	0.0	0.0	20840.546488348853;	%3 onwind
	2	0.0	0.0	3	0.0	0.9985076197246766	0.0	0.0	5689.828737594318;	%3 solar
	2	401333.2651642832	160533.3060657133	3	0.0	2047.6186998177714	0.0	0.0	126610.79544329994;	%4 CCGT
	2	548857.4925934768	219542.9970373907	3	0.0	2800.293329558555	0.0	0.0	57656.27225677763;	%4 OCGT
	2	3807.233667308782	1522.8934669235127	3	0.0	745.2013441590882	0.0	0.0	15872.321081966016;	%4 biomass
	2	0.0	0.0	3	0.0	1.246490381335697	0.0	0.0	20840.546488348853;	%4 onwind
	2	0.0	0.0	3	0.0	0.9772269327335353	0.0	0.0	5689.828737594318;	%4 solar
	2	412221.5688130547	164888.62752522188	3	0.0	2103.1712694543608	0.0	0.0	126610.79544329994;	%5 CCGT
	2	538967.6626472827	215587.06505891308	3	0.0	2749.8350135065443	0.0	0.0	57656.27225677763;	%5 OCGT
	2	3807.183012512747	1522.8732050050987	3	0.0	745.1914293428747	0.0	0.0	15872.321081966016;	%5 biomass
	2	0.0	0.0	3	0.0	1.231336368898412	0.0	0.0	20840.546488348853;	%5 onwind
	2	0.0	0.0	3	0.0	1.0008502728959738	0.0	0.0	5689.828737594318;	%5 solar
	2	404179.20205903857	161671.68082361543	3	0.0	2062.138786015503	0.0	0.0	126610.79544329994;	%6 CCGT
	2	558218.0909823446	223287.23639293783	3	0.0	2848.051484603799	0.0	0.0	57656.27225677763;	%6 OCGT
	2	3807.246143282647	1522.8984573130588	3	0.0	745.2037861191324	0.0	0.0	15872.321081966016;	%6 biomass
	2	0.0	0.0	3	0.0	1.2504845652899788	0.0	0.0	20840.546488348853;	%6 onwind
	2	0.0	0.0	3	0.0	1.0209428451976155	0.0	0.0	5689.828737594318;	%6 solar
	2	383103.3692438105	153241.3476975242	3	0.0	1954.6090267541351	0.0	0.0	126610.79544329994;	%7 CCGT
	2	555831.5921632582	222332.63686530327	3	0.0	2835.8754702207048	0.0	0.0	57656.27225677763;	%7 OCGT
	2	3807.242066237292	1522.896826494917	3	0.0	745.2029881067316	0.0	0.0	15872.321081966016;	%7 biomass
	2	0.0	0.0	3	0.0	1.2641622768738559	0.0	0.0	20840.546488348853;	%7 onwind
	2	0.0	0.0	3	0.0	0.9975052043864983	0.0	0.0	5689.828737594319;	%7 solar
	2	433152.5047495632	173261.00189982526	3	0.0	2209.9617589263426	0.0	0.0	126610.79544329994;	%8 CCGT
	2	568027.7136044257	227211.08544177026	3	0.0	2898.100579614417	0.0	0.0	57656.27225677763;	%8 OCGT
	2	3807.249303562098	1522.899721424839	3	0.0	745.2044046901738	0.0	0.0	15872.321081966016;	%8 biomass
	2	0.0	0.0	3	0.0	1.2778837371516245	0.0	0.0	20840.546488348853;	%8 onwind
	2	0.0	0.0	3	0.0	0.996231944211552	0.0	0.0	5689.828737594318;	%8 solar
	2	415356.5280503558	166142.6112201423	3	0.0	2119.1659594405905	0.0	0.0	126610.79544329994;	%9 CCGT
	2	560430.3205227557	224172.12820910232	3	0.0	2859.33837001406	0.0	0.0	57656.27225677763;	%9 OCGT
	2	3807.0998865437205	1522.8399546174883	3	0.0	745.1751588459034	0.0	0.0	15872.321081966016;	%9 biomass
	2	0.0	0.0	3	0.0	1.2427866089346211	0.0	0.0	20840.546488348853;	%9 onwind
	2	0.0	0.0	3	0.0	1.007842792412104	0.0	0.0	5689.828737594318;	%9 solar
	2	0.0	0.0	3	0.0	0.4834043463273641	0.0	0.0	1012.2443870975252;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.3625532597455231	0.0	0.0	1012.2443870975252;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5266598194714975	0.0	0.0	1012.2443870975252;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.39499486460362315	0.0	0.0	1012.2443870975252;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5495412803214907	0.0	0.0	1012.2443870975252;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.41215596024111806	0.0	0.0	1012.2443870975252;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4589222047765351	0.0	0.0	1012.2443870975252;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.34419165358240134	0.0	0.0	1012.2443870975252;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5459843349556173	0.0	0.0	1012.2443870975252;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.40948825121671295	0.0	0.0	1012.2443870975252;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5252372997847111	0.0	0.0	1012.2443870975252;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.39392797483853337	0.0	0.0	1012.2443870975252;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.499956189134179	0.0	0.0	1012.2443870975252;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.37496714185063423	0.0	0.0	1012.2443870975252;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5074929991075252	0.0	0.0	1012.2443870975252;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.5017809959088319	0.0	0.0	1012.2443870975252;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.37633574693162386	0.0	0.0	1012.2443870975252;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.48862176616803765	0.0	0.0	1012.2443870975252;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.36646632462602824	0.0	0.0	1012.2443870975252;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.47964074099338977	0.0	0.0	1012.2443870975252;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.3597305557450423	0.0	0.0	1012.2443870975252;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4593262645220965	0.0	0.0	1012.2443870975252;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.3444946983915724	0.0	0.0	1012.2443870975252;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5066226042092994	0.0	0.0	1012.2443870975252;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.37996695315697454	0.0	0.0	1012.2443870975252;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4650547961852009	0.0	0.0	1012.2443870975252;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.34879109713890066	0.0	0.0	1012.2443870975252;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46785266124262936	0.0	0.0	1012.2443870975252;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.35088949593197205	0.0	0.0	1012.2443870975252;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5421817544212741	0.0	0.0	1012.2443870975252;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.4066363158159556	0.0	0.0	1012.2443870975252;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4963253151145683	0.0	0.0	1012.2443870975252;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.3722439863359262	0.0	0.0	1012.2443870975252;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5158346290986008	0.0	0.0	1012.2443870975252;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.3868759718239506	0.0	0.0	1012.2443870975252;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5229521137453988	0.0	0.0	1012.2443870975252;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.3922140853090491	0.0	0.0	1012.2443870975252;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46431966545417824	0.0	0.0	1012.2443870975252;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.5404631619353785	0.0	0.0	1012.2443870975252;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.4053473714515339	0.0	0.0	1012.2443870975252;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_012_008
%MPC_012_008	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 12 	Weight: 33
%	Time step: 8

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 33;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	650.9683857388693	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5240.207024572043	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	2231.1216021562896	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	3370.729369041908	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	-686.9103076999029	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	6628.367899531209	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	324.74268476465227	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3024.960115025909	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	3916.8973352804924	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	5397.879887866813	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	166.19740302317595	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	6552.994113166067	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6118.355599787327	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3870.914246564584	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3165.0240161289994	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	3319.955779070805	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8626.532141561302	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5875.102982675016	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3864.5346233611035	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3921.3530936588145	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.316155255054536	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.65050513732342	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.79485388989516	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.384394663553962	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.587471157822325	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.136813875979957	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.09364749695994	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.16371149896221424	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.76543076565658	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.347890849934007	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.03233779726206	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.86920421630329	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.50112202587194	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.46011259454092	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.3444059171575831	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.02202191309294	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.49108150857342	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.314583115515095	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.48388685628513	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.048662975990545	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.42973207412619	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.024954992920854	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	2512.796026311312	1005.1184105245246	3	0.0	491.8371552772189	0.0	0.0	10475.731914097572;	%0 biomass
	2	0.0	0.0	3	0.0	0.8142057769958136	0.0	0.0	12132.550198469726;	%0 offwind
	2	0.0	0.0	3	0.0	0.7952169135400523	0.0	0.0	13754.760682310241;	%0 onwind
	2	0.0	0.0	3	0.0	0.6289424033184126	0.0	0.0	3755.2869668122503;	%0 solar
	2	293135.7921438679	117254.31685754715	3	0.0	1495.5907762442237	0.0	0.0	83563.12499257796;	%1 CCGT
	2	357818.9231500844	143127.56926003375	3	0.0	1825.6067507657365	0.0	0.0	38053.13968947323;	%1 OCGT
	2	2512.6956451243877	1005.0782580497553	3	0.0	491.8175073643351	0.0	0.0	10475.73191409757;	%1 biomass
	2	0.0	0.0	3	0.0	0.8111404617326404	0.0	0.0	13754.760682310241;	%1 onwind
	2	0.0	0.0	3	0.0	0.6580123651069719	0.0	0.0	3755.2869668122503;	%1 solar
	2	256848.14612091423	102739.25844836568	3	0.0	1310.449725106705	0.0	0.0	83563.12499257796;	%10 CCGT
	2	372576.84547225886	149030.73818890358	3	0.0	1900.9022728176474	0.0	0.0	38053.13968947323;	%10 OCGT
	2	2512.741868604201	1005.0967474416803	3	0.0	491.8265548256412	0.0	0.0	10475.73191409757;	%10 biomass
	2	0.0	0.0	3	0.0	0.8148719218663038	0.0	0.0	11977.607560626497;	%10 offwind
	2	0.0	0.0	3	0.0	0.8183904058944492	0.0	0.0	13754.760682310243;	%10 onwind
	2	0.0	0.0	3	0.0	0.6541351743010266	0.0	0.0	3755.2869668122503;	%10 solar
	2	239665.14414051143	95866.05765620458	3	0.0	1222.7813476556705	0.0	0.0	83563.12499257796;	%11 CCGT
	2	355783.06107715884	142313.22443086354	3	0.0	1815.2196993732591	0.0	0.0	38053.13968947323;	%11 OCGT
	2	2512.797440230154	1005.1189760920615	3	0.0	491.83743202782415	0.0	0.0	10475.73191409757;	%11 biomass
	2	0.0	0.0	3	0.0	0.8353692896481932	0.0	0.0	11909.044591400494;	%11 offwind
	2	0.0	0.0	3	0.0	0.7955247398256864	0.0	0.0	13754.760682310241;	%11 onwind
	2	0.0	0.0	3	0.0	0.663148284458847	0.0	0.0	3755.2869668122503;	%11 solar
	2	353865.7198400399	141546.28793601596	3	0.0	1805.4373461226521	0.0	0.0	38053.13968947323;	%12 OCGT
	2	2512.7496828843136	1005.0998731537255	3	0.0	491.8280843382881	0.0	0.0	10475.73191409757;	%12 biomass
	2	0.0	0.0	3	0.0	0.8255264461350277	0.0	0.0	13754.760682310241;	%12 onwind
	2	0.0	0.0	3	0.0	0.6603335390059654	0.0	0.0	3755.2869668122503;	%12 solar
	2	266041.07889292395	106416.43155716959	3	0.0	1357.3524433312446	0.0	0.0	83563.12499257796;	%13 CCGT
	2	370074.83507392963	148029.93402957186	3	0.0	1888.136913642498	0.0	0.0	38053.13968947323;	%13 OCGT
	2	2512.7266173736316	1005.0906469494527	3	0.0	491.82356965622074	0.0	0.0	10475.73191409757;	%13 biomass
	2	0.0	0.0	3	0.0	0.8391648785843179	0.0	0.0	13754.760682310241;	%13 onwind
	2	0.0	0.0	3	0.0	0.6680401915770487	0.0	0.0	3755.2869668122503;	%13 solar
	2	274153.9680640096	109661.58722560384	3	0.0	1398.744735020457	0.0	0.0	83563.12499257796;	%14 CCGT
	2	369809.88433043554	147923.9537321742	3	0.0	1886.7851241348753	0.0	0.0	38053.13968947323;	%14 OCGT
	2	2512.623624310062	1005.0494497240247	3	0.0	491.80341051283256	0.0	0.0	10475.73191409757;	%14 biomass
	2	0.0	0.0	3	0.0	0.8159154417297959	0.0	0.0	12542.334804690077;	%14 offwind
	2	0.0	0.0	3	0.0	0.8078400825715514	0.0	0.0	13754.760682310243;	%14 onwind
	2	0.0	0.0	3	0.0	0.6458055178691516	0.0	0.0	3755.2869668122503;	%14 solar
	2	246214.43562507923	98485.7742500317	3	0.0	1256.196100127955	0.0	0.0	83563.12499257796;	%15 CCGT
	2	345338.59863935085	138135.43945574036	3	0.0	1761.9316257109738	0.0	0.0	38053.13968947323;	%15 OCGT
	2	2512.7552621227587	1005.1021048491035	3	0.0	491.8291763794791	0.0	0.0	10475.73191409757;	%15 biomass
	2	0.0	0.0	3	0.0	0.8260549909603849	0.0	0.0	13089.17743397209;	%15 offwind
	2	0.0	0.0	3	0.0	0.8176980997224144	0.0	0.0	13754.760682310243;	%15 onwind
	2	0.0	0.0	3	0.0	0.638200978928839	0.0	0.0	3755.2869668122503;	%15 solar
	2	257731.07496292752	103092.42998517101	3	0.0	1314.954464096569	0.0	0.0	83563.12499257796;	%16 CCGT
	2	356602.49896151625	142640.9995846065	3	0.0	1819.400504905695	0.0	0.0	38053.13968947323;	%16 OCGT
	2	2512.6947129561486	1005.0778851824595	3	0.0	491.81732490823026	0.0	0.0	10475.731914097572;	%16 biomass
	2	0.0	0.0	3	0.0	0.8163131853591156	0.0	0.0	13754.760682310241;	%16 onwind
	2	0.0	0.0	3	0.0	0.6564464780320448	0.0	0.0	3755.2869668122507;	%16 solar
	2	264981.14148196153	105992.45659278463	3	0.0	1351.9445993977629	0.0	0.0	83563.12499257796;	%17 CCGT
	2	347630.21264463966	139052.08505785585	3	0.0	1773.6235339012226	0.0	0.0	38053.13968947323;	%17 OCGT
	2	2512.7230881864384	1005.0892352745753	3	0.0	491.82287887775266	0.0	0.0	10475.73191409757;	%17 biomass
	2	0.0	0.0	3	0.0	0.8412397611556641	0.0	0.0	13754.760682310243;	%17 onwind
	2	0.0	0.0	3	0.0	0.6726717183646302	0.0	0.0	3755.2869668122507;	%17 solar
	2	277759.9565254137	111103.98261016549	3	0.0	1417.1426353337433	0.0	0.0	83563.12499257796;	%18 CCGT
	2	350335.6819621133	140134.27278484532	3	0.0	1787.4269487862923	0.0	0.0	38053.13968947323;	%18 OCGT
	2	2512.768568238509	1005.1074272954038	3	0.0	491.8317808257015	0.0	0.0	10475.73191409757;	%18 biomass
	2	0.0	0.0	3	0.0	0.8220166444364615	0.0	0.0	13754.760682310241;	%18 onwind
	2	0.0	0.0	3	0.0	0.659981329600072	0.0	0.0	3755.2869668122503;	%18 solar
	2	260264.95393941394	104105.98157576557	3	0.0	1327.8824180582344	0.0	0.0	83563.12499257796;	%19 CCGT
	2	364654.7233062179	145861.88932248717	3	0.0	1860.4832821745808	0.0	0.0	38053.13968947323;	%19 OCGT
	2	2512.7106691611557	1005.0842676644625	3	0.0	491.82044806442667	0.0	0.0	10475.73191409757;	%19 biomass
	2	0.0	0.0	3	0.0	0.8386631834285114	0.0	0.0	13754.760682310241;	%19 onwind
	2	0.0	0.0	3	0.0	0.656012406236478	0.0	0.0	3755.2869668122503;	%19 solar
	2	278638.2037510526	111455.28150042107	3	0.0	1421.6234885257788	0.0	0.0	83563.12499257796;	%2 CCGT
	2	368255.1336611159	147302.05346444636	3	0.0	1878.8527227607951	0.0	0.0	38053.13968947323;	%2 OCGT
	2	2512.69838053003	1005.0793522120122	3	0.0	491.81804277354286	0.0	0.0	10475.73191409757;	%2 biomass
	2	0.0	0.0	3	0.0	0.8264970267957399	0.0	0.0	13754.760682310241;	%2 onwind
	2	0.0	0.0	3	0.0	0.6581196045017396	0.0	0.0	3755.2869668122503;	%2 solar
	2	303861.41132807033	121544.56453122814	3	0.0	1550.3133231023996	0.0	0.0	83563.12499257796;	%3 CCGT
	2	377980.9375902001	151192.37503608002	3	0.0	1928.4741713785718	0.0	0.0	38053.13968947323;	%3 OCGT
	2	2512.624980155999	1005.0499920623996	3	0.0	491.80367589665275	0.0	0.0	10475.73191409757;	%3 biomass
	2	0.0	0.0	3	0.0	0.817552134241991	0.0	0.0	13030.900416554661;	%3 offwind
	2	0.0	0.0	3	0.0	0.8410682630875499	0.0	0.0	13754.760682310243;	%3 onwind
	2	0.0	0.0	3	0.0	0.6590150290182866	0.0	0.0	3755.2869668122503;	%3 solar
	2	264879.9550084269	105951.98200337077	3	0.0	1351.428341879729	0.0	0.0	83563.12499257796;	%4 CCGT
	2	362245.9451116947	144898.37804467787	3	0.0	1848.1935975086462	0.0	0.0	38053.13968947323;	%4 OCGT
	2	2512.774220423796	1005.1096881695183	3	0.0	491.8328871449982	0.0	0.0	10475.73191409757;	%4 biomass
	2	0.0	0.0	3	0.0	0.82268365168156	0.0	0.0	13754.760682310241;	%4 onwind
	2	0.0	0.0	3	0.0	0.6449697756041333	0.0	0.0	3755.2869668122503;	%4 solar
	2	272066.23541661614	108826.49416664644	3	0.0	1388.093037839878	0.0	0.0	83563.12499257796;	%5 CCGT
	2	355718.6573472066	142287.46293888264	3	0.0	1814.8911089143194	0.0	0.0	38053.13968947323;	%5 OCGT
	2	2512.740788258413	1005.0963153033653	3	0.0	491.82634336629735	0.0	0.0	10475.73191409757;	%5 biomass
	2	0.0	0.0	3	0.0	0.8126820034729519	0.0	0.0	13754.760682310241;	%5 onwind
	2	0.0	0.0	3	0.0	0.6605611801113428	0.0	0.0	3755.2869668122503;	%5 solar
	2	266758.27335896547	106703.30934358618	3	0.0	1361.0115987702318	0.0	0.0	83563.12499257796;	%6 CCGT
	2	368423.94004834746	147369.57601933897	3	0.0	1879.7139798385074	0.0	0.0	38053.13968947323;	%6 OCGT
	2	2512.782454566547	1005.1129818266188	3	0.0	491.8344988386273	0.0	0.0	10475.73191409757;	%6 biomass
	2	0.0	0.0	3	0.0	0.825319813091386	0.0	0.0	13754.760682310241;	%6 onwind
	2	0.0	0.0	3	0.0	0.6738222778304263	0.0	0.0	3755.2869668122503;	%6 solar
	2	252848.22370091494	101139.28948036597	3	0.0	1290.0419576577292	0.0	0.0	83563.12499257796;	%7 CCGT
	2	366848.8508277504	146739.54033110014	3	0.0	1871.6778103456652	0.0	0.0	38053.13968947323;	%7 OCGT
	2	2512.7797637166127	1005.1119054866451	3	0.0	491.8339721504429	0.0	0.0	10475.73191409757;	%7 biomass
	2	0.0	0.0	3	0.0	0.8343471027367448	0.0	0.0	13754.760682310241;	%7 onwind
	2	0.0	0.0	3	0.0	0.6583534348950889	0.0	0.0	3755.2869668122507;	%7 solar
	2	285880.6531347117	114352.26125388467	3	0.0	1458.574760891386	0.0	0.0	83563.12499257796;	%8 CCGT
	2	374898.29097892094	149959.31639156837	3	0.0	1912.746382545515	0.0	0.0	38053.13968947323;	%8 OCGT
	2	2512.7845403509846	1005.1138161403939	3	0.0	491.83490709551467	0.0	0.0	10475.73191409757;	%8 biomass
	2	0.0	0.0	3	0.0	0.8434032665200721	0.0	0.0	13754.760682310241;	%8 onwind
	2	0.0	0.0	3	0.0	0.6575130831796243	0.0	0.0	3755.2869668122503;	%8 solar
	2	274135.30851323484	109654.12340529393	3	0.0	1398.6495332307898	0.0	0.0	83563.12499257796;	%9 CCGT
	2	369884.0115450188	147953.60461800752	3	0.0	1887.1633242092796	0.0	0.0	38053.13968947323;	%9 OCGT
	2	2512.6859251188553	1005.0743700475423	3	0.0	491.8156048382963	0.0	0.0	10475.73191409757;	%9 biomass
	2	0.0	0.0	3	0.0	0.82023916189685	0.0	0.0	13754.760682310241;	%9 onwind
	2	0.0	0.0	3	0.0	0.6651762429919887	0.0	0.0	3755.2869668122503;	%9 solar
	2	0.0	0.0	3	0.0	0.31904686857606035	0.0	0.0	668.0812954843666;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.23928515143204526	0.0	0.0	668.0812954843666;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34759548085118835	0.0	0.0	668.0812954843666;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.2606966106383913	0.0	0.0	668.0812954843666;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3626972450121839	0.0	0.0	668.0812954843666;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.27202293375913794	0.0	0.0	668.0812954843666;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3028886551525132	0.0	0.0	668.0812954843666;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.22716649136438488	0.0	0.0	668.0812954843666;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3603496610707074	0.0	0.0	668.0812954843666;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.2702622458030306	0.0	0.0	668.0812954843666;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34665661785790935	0.0	0.0	668.0812954843666;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.259992463393432	0.0	0.0	668.0812954843666;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32997108482855814	0.0	0.0	668.0812954843666;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.2474783136214186	0.0	0.0	668.0812954843666;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33494537941096664	0.0	0.0	668.0812954843666;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.33117545729982906	0.0	0.0	668.0812954843666;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.24838159297487178	0.0	0.0	668.0812954843666;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32249036567090483	0.0	0.0	668.0812954843666;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.24186777425317862	0.0	0.0	668.0812954843666;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31656288905563723	0.0	0.0	668.0812954843666;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.23742216679172792	0.0	0.0	668.0812954843666;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3031553345845837	0.0	0.0	668.0812954843666;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.2273665009384378	0.0	0.0	668.0812954843666;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3343709187781376	0.0	0.0	668.0812954843666;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.2507781890836032	0.0	0.0	668.0812954843666;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3069361654822326	0.0	0.0	668.0812954843666;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.23020212411167446	0.0	0.0	668.0812954843666;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30878275642013536	0.0	0.0	668.0812954843666;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.23158706731510154	0.0	0.0	668.0812954843666;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.35783995791804096	0.0	0.0	668.0812954843666;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.26837996843853074	0.0	0.0	668.0812954843666;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32757470797561505	0.0	0.0	668.0812954843666;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.2456810309817113	0.0	0.0	668.0812954843666;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34045085520507656	0.0	0.0	668.0812954843666;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.2553381414038074	0.0	0.0	668.0812954843666;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3451483950719632	0.0	0.0	668.0812954843666;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.2588612963039724	0.0	0.0	668.0812954843666;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3064509791997576	0.0	0.0	668.0812954843666;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.35670568687734977	0.0	0.0	668.0812954843666;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.26752926515801234	0.0	0.0	668.0812954843666;	%DE1 76 PHS (pump mode)
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

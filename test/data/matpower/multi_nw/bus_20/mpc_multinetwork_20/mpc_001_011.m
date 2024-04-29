function mpc = mpc_001_011
%MPC_001_011	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 1 	Weight: 31
%	Time step: 11

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 31;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	913.8772491449608	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6373.373907012776	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1997.2796430347164	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	1044.7052608870067	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	697.1343741416056	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	7909.872521988441	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	390.0426398432281	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	3448.798045886102	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2618.126659766334	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	3650.355533387635	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	3004.551699437215	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	5777.509684387451	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6712.941142258053	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3879.874484374046	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2485.763435730243	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-1009.3736652194747	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	7468.218681544805	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	6079.284226668193	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4412.135876925713	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	5315.770068127738	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.9235314485377835	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.489101349095055	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.804264444707526	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.555969727867653	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.057669704427215	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.1928441651946327	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.604526172426555	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.86778355899011	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.167752197205288	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.346566440360999	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.410133446194736	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.334905369751956	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.285136844017304	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.832748201345535	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.71840770763116	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.142690904554841	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.7797538235461925	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.16388526410062	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.471247587821565	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.260583011640665	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.142307661356698	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.955549289178118	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.2728988083037445	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.148142407337117	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.4640866046823415	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.597785183706005	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.756284451928566	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.073419682328495	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.149799940176135	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.777994117343548	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.933978670500957	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.55521407719683	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.253811242266808	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.688809128539702	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.628718616433116	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.04751809755796	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.93086721785794	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.737502361060807	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.348793472763545	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.687040038929504	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.780518349588642	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.380289705049215	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.8511454060962915	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.571923850374857	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.073106976934614	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.606126815115075	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	2360.50535805002	944.202143220008	3	0.0	462.0288428361754	0.0	0.0	9840.83907081893;	%0 biomass
	2	0.0	0.0	3	0.0	0.7648599723294006	0.0	0.0	11397.244125835197;	%0 offwind
	2	0.0	0.0	3	0.0	0.7470219490830795	0.0	0.0	12921.138822776289;	%0 onwind
	2	0.0	0.0	3	0.0	0.5908246819051755	0.0	0.0	3527.6938173084773;	%0 solar
	2	275369.9865593911	110147.99462375641	3	0.0	1404.948911017301	0.0	0.0	78498.69317484596;	%1 CCGT
	2	336132.927807655	134453.171123062	3	0.0	1714.9639173859948	0.0	0.0	35746.888799202134;	%1 OCGT
	2	2360.411060571395	944.1644242285579	3	0.0	462.01038570589054	0.0	0.0	9840.839070818929;	%1 biomass
	2	0.0	0.0	3	0.0	0.7619804337488439	0.0	0.0	12921.138822776289;	%1 onwind
	2	0.0	0.0	3	0.0	0.6181328278277615	0.0	0.0	3527.6938173084773;	%1 solar
	2	241281.5918105558	96512.63672422232	3	0.0	1231.0285296456927	0.0	0.0	78498.69317484596;	%10 CCGT
	2	349996.4305951523	139998.5722380609	3	0.0	1785.6960744650628	0.0	0.0	35746.888799202134;	%10 OCGT
	2	2360.454482628189	944.1817930512756	3	0.0	462.0188848362084	0.0	0.0	9840.839070818929;	%10 biomass
	2	0.0	0.0	3	0.0	0.7654857447834974	0.0	0.0	11251.691950891556;	%10 offwind
	2	0.0	0.0	3	0.0	0.7687909873553916	0.0	0.0	12921.13882277629;	%10 onwind
	2	0.0	0.0	3	0.0	0.6144906182827825	0.0	0.0	3527.6938173084773;	%10 solar
	2	225139.98388957136	90055.99355582854	3	0.0	1148.6733871916906	0.0	0.0	78498.69317484596;	%11 CCGT
	2	334220.45131490676	133688.1805259627	3	0.0	1705.206384259728	0.0	0.0	35746.888799202134;	%11 OCGT
	2	2360.5066862768113	944.2026745107244	3	0.0	462.02910281401665	0.0	0.0	9840.839070818929;	%11 biomass
	2	0.0	0.0	3	0.0	0.7847408478513331	0.0	0.0	11187.284313133796;	%11 offwind
	2	0.0	0.0	3	0.0	0.7473111192301902	0.0	0.0	12921.138822776289;	%11 onwind
	2	0.0	0.0	3	0.0	0.6229574793401289	0.0	0.0	3527.6938173084773;	%11 solar
	2	332419.31257700716	132967.72503080286	3	0.0	1696.0169009030976	0.0	0.0	35746.888799202134;	%12 OCGT
	2	2360.4618233155675	944.1847293262271	3	0.0	462.0203216511191	0.0	0.0	9840.839070818929;	%12 biomass
	2	0.0	0.0	3	0.0	0.7754945403086623	0.0	0.0	12921.138822776289;	%12 onwind
	2	0.0	0.0	3	0.0	0.6203133245207554	0.0	0.0	3527.6938173084773;	%12 solar
	2	249917.37714183767	99966.95085673506	3	0.0	1275.0886588869266	0.0	0.0	78498.69317484596;	%13 CCGT
	2	347646.0571906612	139058.42287626446	3	0.0	1773.7043734217407	0.0	0.0	35746.888799202134;	%13 OCGT
	2	2360.440155714624	944.1760622858495	3	0.0	462.01608058614676	0.0	0.0	9840.839070818929;	%13 biomass
	2	0.0	0.0	3	0.0	0.7883064010943591	0.0	0.0	12921.138822776289;	%13 onwind
	2	0.0	0.0	3	0.0	0.6275529072390459	0.0	0.0	3527.6938173084773;	%13 solar
	2	257538.57606013023	103015.43042405209	3	0.0	1313.9723268373991	0.0	0.0	78498.69317484596;	%14 CCGT
	2	347397.1640679849	138958.86562719397	3	0.0	1772.4345105509433	0.0	0.0	35746.888799202134;	%14 OCGT
	2	2360.3434046549064	944.1373618619625	3	0.0	461.99714320902456	0.0	0.0	9840.839070818929;	%14 biomass
	2	0.0	0.0	3	0.0	0.7664660210188993	0.0	0.0	11782.193301375528;	%14 offwind
	2	0.0	0.0	3	0.0	0.758880077567215	0.0	0.0	12921.13882277629;	%14 onwind
	2	0.0	0.0	3	0.0	0.6066657895134454	0.0	0.0	3527.6938173084773;	%14 solar
	2	231292.34861749868	92516.93944699947	3	0.0	1180.0630031505034	0.0	0.0	78498.69317484596;	%15 CCGT
	2	324408.9866006023	129763.59464024093	3	0.0	1655.1478908193997	0.0	0.0	35746.888799202134;	%15 OCGT
	2	2360.467064418349	944.1868257673397	3	0.0	462.02134750799553	0.0	0.0	9840.839070818929;	%15 biomass
	2	0.0	0.0	3	0.0	0.775991052114301	0.0	0.0	12295.893953125298;	%15 offwind
	2	0.0	0.0	3	0.0	0.7681406391331771	0.0	0.0	12921.13882277629;	%15 onwind
	2	0.0	0.0	3	0.0	0.5995221317210306	0.0	0.0	3527.6938173084773;	%15 solar
	2	242111.0098136592	96844.40392546367	3	0.0	1235.2602541513222	0.0	0.0	78498.69317484596;	%16 CCGT
	2	334990.226297182	133996.09051887278	3	0.0	1709.1338076386833	0.0	0.0	35746.888799202134;	%16 OCGT
	2	2360.4101848982	944.16407395928	3	0.0	462.0102143077314	0.0	0.0	9840.83907081893;	%16 biomass
	2	0.0	0.0	3	0.0	0.7668396589737146	0.0	0.0	12921.138822776289;	%16 onwind
	2	0.0	0.0	3	0.0	0.6166618429997998	0.0	0.0	3527.6938173084777;	%16 solar
	2	248921.67836184267	99568.67134473707	3	0.0	1270.0085630706258	0.0	0.0	78498.69317484596;	%17 CCGT
	2	326561.7149086009	130624.68596344035	3	0.0	1666.1311985132695	0.0	0.0	35746.888799202134;	%17 OCGT
	2	2360.4368404175634	944.1747361670253	3	0.0	462.0154316730404	0.0	0.0	9840.839070818929;	%17 biomass
	2	0.0	0.0	3	0.0	0.790255533206836	0.0	0.0	12921.13882277629;	%17 onwind
	2	0.0	0.0	3	0.0	0.6319037354334405	0.0	0.0	3527.6938173084777;	%17 solar
	2	260926.01976629772	104370.4079065191	3	0.0	1331.2552028892742	0.0	0.0	78498.69317484596;	%18 CCGT
	2	329103.2163886519	131641.28655546077	3	0.0	1679.0980427992442	0.0	0.0	35746.888799202134;	%18 OCGT
	2	2360.4795641028422	944.1918256411368	3	0.0	462.0237941089924	0.0	0.0	9840.839070818929;	%18 biomass
	2	0.0	0.0	3	0.0	0.7721974538645547	0.0	0.0	12921.138822776289;	%18 onwind
	2	0.0	0.0	3	0.0	0.6199824611394615	0.0	0.0	3527.6938173084773;	%18 solar
	2	244491.32036732824	97796.5281469313	3	0.0	1247.4046957516746	0.0	0.0	78498.69317484596;	%19 CCGT
	2	342554.437045235	137021.77481809398	3	0.0	1747.7267196185458	0.0	0.0	35746.888799202134;	%19 OCGT
	2	2360.42517406048	944.1700696241919	3	0.0	462.01314818173415	0.0	0.0	9840.839070818929;	%19 biomass
	2	0.0	0.0	3	0.0	0.7878351117055713	0.0	0.0	12921.138822776289;	%19 onwind
	2	0.0	0.0	3	0.0	0.6162540785857824	0.0	0.0	3527.6938173084773;	%19 solar
	2	261751.03988735247	104700.415954941	3	0.0	1335.464489221186	0.0	0.0	78498.69317484596;	%2 CCGT
	2	345936.64071195736	138374.65628478295	3	0.0	1764.9828607752925	0.0	0.0	35746.888799202134;	%2 OCGT
	2	2360.413630194877	944.1654520779508	3	0.0	462.0108886660554	0.0	0.0	9840.839070818929;	%2 biomass
	2	0.0	0.0	3	0.0	0.7764062978990285	0.0	0.0	12921.138822776289;	%2 onwind
	2	0.0	0.0	3	0.0	0.6182335678652705	0.0	0.0	3527.6938173084773;	%2 solar
	2	285445.5682172782	114178.2272869113	3	0.0	1456.3549398840726	0.0	0.0	78498.69317484596;	%3 CCGT
	2	355073.0019786728	142029.20079146913	3	0.0	1811.5969488707794	0.0	0.0	35746.888799202134;	%3 OCGT
	2	2360.3446783283625	944.1378713313451	3	0.0	461.9973925089769	0.0	0.0	9840.839070818929;	%3 biomass
	2	0.0	0.0	3	0.0	0.7680035200455068	0.0	0.0	12241.148876157407;	%3 offwind
	2	0.0	0.0	3	0.0	0.7900944289610317	0.0	0.0	12921.13882277629;	%3 onwind
	2	0.0	0.0	3	0.0	0.6190747242292995	0.0	0.0	3527.6938173084773;	%3 solar
	2	248826.62440185557	99530.64976074224	3	0.0	1269.5235938870183	0.0	0.0	78498.69317484596;	%4 CCGT
	2	340291.64540795557	136116.65816318223	3	0.0	1736.181864326304	0.0	0.0	35746.888799202134;	%4 OCGT
	2	2360.4848737314446	944.1939494925779	3	0.0	462.0248333786347	0.0	0.0	9840.839070818929;	%4 biomass
	2	0.0	0.0	3	0.0	0.7728240364281321	0.0	0.0	12921.138822776289;	%4 onwind
	2	0.0	0.0	3	0.0	0.605880698294792	0.0	0.0	3527.6938173084773;	%4 solar
	2	255577.37266409394	102230.94906563757	3	0.0	1303.9661870617037	0.0	0.0	78498.69317484596;	%5 CCGT
	2	334159.9508413153	133663.98033652612	3	0.0	1704.8977083740574	0.0	0.0	35746.888799202134;	%5 OCGT
	2	2360.453467757903	944.1813871031612	3	0.0	462.0186861925823	0.0	0.0	9840.839070818929;	%5 biomass
	2	0.0	0.0	3	0.0	0.7634285487170155	0.0	0.0	12921.138822776289;	%5 onwind
	2	0.0	0.0	3	0.0	0.6205271691955038	0.0	0.0	3527.6938173084773;	%5 solar
	2	250591.10527660392	100236.44211064157	3	0.0	1278.5260473296119	0.0	0.0	78498.69317484596;	%6 CCGT
	2	346095.21640905365	138438.08656362147	3	0.0	1765.7919204543555	0.0	0.0	35746.888799202134;	%6 OCGT
	2	2360.4926088352413	944.1970435340966	3	0.0	462.02634739386207	0.0	0.0	9840.839070818929;	%6 biomass
	2	0.0	0.0	3	0.0	0.7753004304797869	0.0	0.0	12921.138822776289;	%6 onwind
	2	0.0	0.0	3	0.0	0.6329845640225217	0.0	0.0	3527.6938173084773;	%6 solar
	2	237524.0889311625	95009.635572465	3	0.0	1211.8575965875639	0.0	0.0	78498.69317484596;	%7 CCGT
	2	344615.5871412201	137846.23485648804	3	0.0	1758.2427915368369	0.0	0.0	35746.888799202134;	%7 OCGT
	2	2360.490081067121	944.1960324268484	3	0.0	462.0258526261736	0.0	0.0	9840.839070818929;	%7 biomass
	2	0.0	0.0	3	0.0	0.7837806116617907	0.0	0.0	12921.138822776289;	%7 onwind
	2	0.0	0.0	3	0.0	0.6184532267196289	0.0	0.0	3527.6938173084777;	%7 solar
	2	268554.55294472916	107421.82117789166	3	0.0	1370.1762905343323	0.0	0.0	78498.69317484596;	%8 CCGT
	2	352177.1824347439	140870.87297389758	3	0.0	1796.8223593609384	0.0	0.0	35746.888799202134;	%8 OCGT
	2	2360.494568208501	944.1978272834002	3	0.0	462.02673090790773	0.0	0.0	9840.839070818929;	%8 biomass
	2	0.0	0.0	3	0.0	0.7922879170340071	0.0	0.0	12921.138822776289;	%8 onwind
	2	0.0	0.0	3	0.0	0.6176638054111623	0.0	0.0	3527.6938173084773;	%8 solar
	2	257521.04739122058	103008.41895648824	3	0.0	1313.8828948531661	0.0	0.0	78498.69317484596;	%9 CCGT
	2	347466.79872410855	138986.71948964344	3	0.0	1772.789789408717	0.0	0.0	35746.888799202134;	%9 OCGT
	2	2360.401929657107	944.1607718628427	3	0.0	462.0085984844601	0.0	0.0	9840.839070818929;	%9 biomass
	2	0.0	0.0	3	0.0	0.7705276975394652	0.0	0.0	12921.138822776289;	%9 onwind
	2	0.0	0.0	3	0.0	0.6248625312955045	0.0	0.0	3527.6938173084773;	%9 solar
	2	0.0	0.0	3	0.0	0.29971069472296574	0.0	0.0	627.5915200004656;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.2247830210422243	0.0	0.0	627.5915200004656;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.32652908807232844	0.0	0.0	627.5915200004656;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.24489681605424635	0.0	0.0	627.5915200004656;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34071559379932426	0.0	0.0	627.5915200004656;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.2555366953494932	0.0	0.0	627.5915200004656;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28453176696145177	0.0	0.0	627.5915200004656;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.21339882522108883	0.0	0.0	627.5915200004656;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3385102876724827	0.0	0.0	627.5915200004656;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.25388271575436205	0.0	0.0	627.5915200004656;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3256471258665209	0.0	0.0	627.5915200004656;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.24423534439989067	0.0	0.0	627.5915200004656;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30997283726319097	0.0	0.0	627.5915200004656;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.23247962794739324	0.0	0.0	627.5915200004656;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3146456594466656	0.0	0.0	627.5915200004656;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.3111042174634758	0.0	0.0	627.5915200004656;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.23332816309760684	0.0	0.0	627.5915200004656;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30294549502418333	0.0	0.0	627.5915200004656;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.2272091212681375	0.0	0.0	627.5915200004656;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.29737725941590165	0.0	0.0	627.5915200004656;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.22303294456192624	0.0	0.0	627.5915200004656;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2847822840036998	0.0	0.0	627.5915200004656;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.21358671300277488	0.0	0.0	627.5915200004656;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3141060146097656	0.0	0.0	627.5915200004656;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.2355795109573242	0.0	0.0	627.5915200004656;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28833397363482455	0.0	0.0	627.5915200004656;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.21625048022611842	0.0	0.0	627.5915200004656;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2900686499704302	0.0	0.0	627.5915200004656;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.21755148747782266	0.0	0.0	627.5915200004656;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33615268774118995	0.0	0.0	627.5915200004656;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.25211451580589245	0.0	0.0	627.5915200004656;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30772169537103233	0.0	0.0	627.5915200004656;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.23079127152827425	0.0	0.0	627.5915200004656;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3198174700411325	0.0	0.0	627.5915200004656;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.23986310253084936	0.0	0.0	627.5915200004656;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3242303105221473	0.0	0.0	627.5915200004656;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.24317273289161045	0.0	0.0	627.5915200004656;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28787819258159053	0.0	0.0	627.5915200004656;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.33508716039993464	0.0	0.0	627.5915200004656;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.251315370299951	0.0	0.0	627.5915200004656;	%DE1 76 PHS (pump mode)
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

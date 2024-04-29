function mpc = mpc_003_014
%MPC_003_014	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 3 	Weight: 13
%	Time step: 14

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 13;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1192.0109194785473	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4905.553621286485	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1743.451305253147	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	536.601230573889	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	4516.15673183592	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	6310.295674576133	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	265.6026658278136	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2558.2860676511077	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2572.273322909648	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	1579.293064180582	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1460.5030444747824	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	3182.016080657476	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5322.006216066827	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3747.5146441417323	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	4094.774091258907	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	2687.0973222008693	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3068.3381819475976	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4597.16608429141	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5622.788860444541	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3797.21605350359	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.580212140118102	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.097492357264265	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.606480953845605	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.524083691096031	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.371420047016351	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.148135760829028	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.709643250685573	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.5668852206966415	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.342856237453846	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.876437882359856	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.529213899009432	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.22181101162299	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.55038352654576	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.993641308034643	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.278257600570908	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.758550183432067	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.486666244457403	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.610511033162631	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.152821813785176	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.277163060398136	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6348319002947078	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.371155259725542	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.3855187208433148	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.932241073530584	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6879091738544485	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.924172333587222	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.0967912169096021	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.72617350633881	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6761825187736334	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.35948071083162	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8290658934271865	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.727727434916952	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.258609421405872	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.7354382985674195	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.143287177075731	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.4708475565351848	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.86044440361994	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.7090186295699756	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.36937856228	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.6632182877220794	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.1408759254818	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.1983536901172693	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.530440525100154	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5250542961809014	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.727090945787518	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.374214244871122	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	989.8893436983955	395.95573747935816	3	0.0	193.7540308667832	0.0	0.0	4126.803481311164;	%0 biomass
	2	0.0	0.0	3	0.0	0.32074773033168413	0.0	0.0	4779.489472124437;	%0 offwind
	2	0.0	0.0	3	0.0	0.31326726897032364	0.0	0.0	5418.542086970701;	%0 onwind
	2	0.0	0.0	3	0.0	0.24776518918604132	0.0	0.0	1479.3554717745228;	%0 solar
	2	115477.73629909947	46191.09451963979	3	0.0	589.1721239749972	0.0	0.0	32918.80681525798;	%1 CCGT
	2	140958.9697257908	56383.587890316325	3	0.0	719.1784169683204	0.0	0.0	14990.630786762184;	%1 OCGT
	2	989.8497995944558	395.9399198377824	3	0.0	193.7462907798896	0.0	0.0	4126.803481311164;	%1 biomass
	2	0.0	0.0	3	0.0	0.31954018189467653	0.0	0.0	5418.542086970701;	%1 onwind
	2	0.0	0.0	3	0.0	0.2592169923148677	0.0	0.0	1479.3554717745228;	%1 solar
	2	101182.60301732985	40473.041206931935	3	0.0	516.2377704965808	0.0	0.0	32918.80681525798;	%10 CCGT
	2	146772.6967011929	58709.07868047716	3	0.0	748.8402892918004	0.0	0.0	14990.630786762184;	%10 OCGT
	2	989.8680088440792	395.94720353763165	3	0.0	193.7498549313132	0.0	0.0	4126.803481311164;	%10 biomass
	2	0.0	0.0	3	0.0	0.3210101510382409	0.0	0.0	4718.451463277104;	%10 offwind
	2	0.0	0.0	3	0.0	0.32239622050387395	0.0	0.0	5418.542086970702;	%10 onwind
	2	0.0	0.0	3	0.0	0.2576896141185862	0.0	0.0	1479.3554717745228;	%10 solar
	2	94413.54163111057	37765.416652444226	3	0.0	481.70174301587025	0.0	0.0	32918.80681525798;	%11 CCGT
	2	140156.96345463832	56062.78538185533	3	0.0	715.0865482379505	0.0	0.0	14990.630786762184;	%11 OCGT
	2	989.8899006967273	395.9559602786909	3	0.0	193.7541398897489	0.0	0.0	4126.803481311164;	%11 biomass
	2	0.0	0.0	3	0.0	0.3290848716795913	0.0	0.0	4691.441808733528;	%11 offwind
	2	0.0	0.0	3	0.0	0.31338853387072496	0.0	0.0	5418.542086970701;	%11 onwind
	2	0.0	0.0	3	0.0	0.26124023327166696	0.0	0.0	1479.3554717745228;	%11 solar
	2	139401.64720971268	55760.65888388507	3	0.0	711.2328939271055	0.0	0.0	14990.630786762184;	%12 OCGT
	2	989.8710871968508	395.9484348787404	3	0.0	193.75045746659833	0.0	0.0	4126.803481311164;	%12 biomass
	2	0.0	0.0	3	0.0	0.32520738787137454	0.0	0.0	5418.542086970701;	%12 onwind
	2	0.0	0.0	3	0.0	0.26013139415386516	0.0	0.0	1479.3554717745228;	%12 solar
	2	104804.06138206096	41921.62455282438	3	0.0	534.714598888066	0.0	0.0	32918.80681525798;	%13 CCGT
	2	145787.056241245	58314.822496498	3	0.0	743.8115114349234	0.0	0.0	14990.630786762184;	%13 OCGT
	2	989.8620007835519	395.9448003134207	3	0.0	193.7486789554809	0.0	0.0	4126.803481311164;	%13 biomass
	2	0.0	0.0	3	0.0	0.3305801036847313	0.0	0.0	5418.542086970701;	%13 onwind
	2	0.0	0.0	3	0.0	0.2631673481970192	0.0	0.0	1479.3554717745228;	%13 solar
	2	108000.0480252159	43200.019210086364	3	0.0	551.0206531898771	0.0	0.0	32918.80681525798;	%14 CCGT
	2	145682.68170592916	58273.07268237166	3	0.0	743.2789882955569	0.0	0.0	14990.630786762184;	%14 OCGT
	2	989.8214277585092	395.92857110340367	3	0.0	193.74073747475222	0.0	0.0	4126.803481311164;	%14 biomass
	2	0.0	0.0	3	0.0	0.3214212346208287	0.0	0.0	4940.919771544576;	%14 offwind
	2	0.0	0.0	3	0.0	0.3182400325281869	0.0	0.0	5418.542086970702;	%14 onwind
	2	0.0	0.0	3	0.0	0.25440823431209003	0.0	0.0	1479.3554717745228;	%14 solar
	2	96993.56554927364	38797.42621970946	3	0.0	494.8651303534369	0.0	0.0	32918.80681525798;	%15 CCGT
	2	136042.47825186548	54416.991300746195	3	0.0	694.0942767952321	0.0	0.0	14990.630786762184;	%15 OCGT
	2	989.8732850786624	395.949314031465	3	0.0	193.7508876646433	0.0	0.0	4126.803481311164;	%15 biomass
	2	0.0	0.0	3	0.0	0.3254156024995456	0.0	0.0	5156.342625504157;	%15 offwind
	2	0.0	0.0	3	0.0	0.322123493830042	0.0	0.0	5418.542086970702;	%15 onwind
	2	0.0	0.0	3	0.0	0.2514125068507548	0.0	0.0	1479.3554717745228;	%15 solar
	2	101530.42347024418	40612.16938809767	3	0.0	518.012364644103	0.0	0.0	32918.80681525798;	%16 CCGT
	2	140479.77231817308	56191.90892726923	3	0.0	716.7335322355768	0.0	0.0	14990.630786762184;	%16 OCGT
	2	989.8494323766646	395.9397729506658	3	0.0	193.7462189032422	0.0	0.0	4126.803481311164;	%16 biomass
	2	0.0	0.0	3	0.0	0.3215779215051061	0.0	0.0	5418.542086970701;	%16 onwind
	2	0.0	0.0	3	0.0	0.25860012770959345	0.0	0.0	1479.355471774523;	%16 solar
	2	104386.51028077274	41754.6041123091	3	0.0	532.5842361263915	0.0	0.0	32918.80681525798;	%17 CCGT
	2	136945.23528425198	54778.09411370079	3	0.0	698.7001800216937	0.0	0.0	14990.630786762184;	%17 OCGT
	2	989.8606104976878	395.9442441990751	3	0.0	193.74840683062985	0.0	0.0	4126.803481311164;	%17 biomass
	2	0.0	0.0	3	0.0	0.3313974816673828	0.0	0.0	5418.542086970702;	%17 onwind
	2	0.0	0.0	3	0.0	0.26499188905273313	0.0	0.0	1479.355471774523;	%17 solar
	2	109420.58893425389	43768.23557370156	3	0.0	558.2683108890504	0.0	0.0	32918.80681525798;	%18 CCGT
	2	138011.02622749918	55204.410490999675	3	0.0	704.1378889158121	0.0	0.0	14990.630786762184;	%18 OCGT
	2	989.8785268818369	395.95141075273483	3	0.0	193.7519136586097	0.0	0.0	4126.803481311164;	%18 biomass
	2	0.0	0.0	3	0.0	0.3238247387173939	0.0	0.0	5418.542086970701;	%18 onwind
	2	0.0	0.0	3	0.0	0.25999264499396774	0.0	0.0	1479.3554717745228;	%18 solar
	2	102528.618218557	41011.4472874228	3	0.0	523.1051949926377	0.0	0.0	32918.80681525798;	%19 CCGT
	2	143651.86069638887	57460.74427855555	3	0.0	732.9176566142288	0.0	0.0	14990.630786762184;	%19 OCGT
	2	989.8557181543947	395.94228726175794	3	0.0	193.74744923750143	0.0	0.0	4126.803481311164;	%19 biomass
	2	0.0	0.0	3	0.0	0.33038246619911055	0.0	0.0	5418.542086970701;	%19 onwind
	2	0.0	0.0	3	0.0	0.25842912972952164	0.0	0.0	1479.3554717745228;	%19 solar
	2	109766.56511405105	43906.62604562042	3	0.0	560.0334954798523	0.0	0.0	32918.80681525798;	%2 CCGT
	2	145070.2041695305	58028.08166781221	3	0.0	740.1541029057678	0.0	0.0	14990.630786762184;	%2 OCGT
	2	989.8508771784967	395.9403508713987	3	0.0	193.7465016986684	0.0	0.0	4126.803481311164;	%2 biomass
	2	0.0	0.0	3	0.0	0.3255897378286248	0.0	0.0	5418.542086970701;	%2 onwind
	2	0.0	0.0	3	0.0	0.25925923813704893	0.0	0.0	1479.3554717745228;	%2 solar
	2	119702.98022014892	47881.19208805957	3	0.0	610.7294909191272	0.0	0.0	32918.80681525798;	%3 CCGT
	2	148901.5814749273	59560.63258997092	3	0.0	759.7019463006495	0.0	0.0	14990.630786762184;	%3 OCGT
	2	989.8219618796359	395.9287847518544	3	0.0	193.7408420198935	0.0	0.0	4126.803481311164;	%3 biomass
	2	0.0	0.0	3	0.0	0.322065992277148	0.0	0.0	5133.385012582139;	%3 offwind
	2	0.0	0.0	3	0.0	0.33132992182236815	0.0	0.0	5418.542086970702;	%3 onwind
	2	0.0	0.0	3	0.0	0.25961198112841594	0.0	0.0	1479.3554717745228;	%3 solar
	2	104346.64894271363	41738.659577085455	3	0.0	532.3808619526205	0.0	0.0	32918.80681525798;	%4 CCGT
	2	142702.94807430395	57081.17922972159	3	0.0	728.0762656852243	0.0	0.0	14990.630786762184;	%4 OCGT
	2	989.8807535002833	395.9523014001133	3	0.0	193.75234948136293	0.0	0.0	4126.803481311164;	%4 biomass
	2	0.0	0.0	3	0.0	0.3240874991472812	0.0	0.0	5418.542086970701;	%4 onwind
	2	0.0	0.0	3	0.0	0.2540790025107192	0.0	0.0	1479.3554717745228;	%4 solar
	2	107177.60789139423	42871.04315655769	3	0.0	546.8245300581337	0.0	0.0	32918.80681525798;	%5 CCGT
	2	140131.59228829353	56052.6369153174	3	0.0	714.9571035117016	0.0	0.0	14990.630786762184;	%5 OCGT
	2	989.8675832533143	395.9470333013257	3	0.0	193.74977162914743	0.0	0.0	4126.803481311164;	%5 biomass
	2	0.0	0.0	3	0.0	0.3201474559135871	0.0	0.0	5418.542086970701;	%5 onwind
	2	0.0	0.0	3	0.0	0.2602210709529532	0.0	0.0	1479.3554717745228;	%5 solar
	2	105086.59253535002	42034.63701414001	3	0.0	536.1560843640308	0.0	0.0	32918.80681525798;	%6 CCGT
	2	145136.7036554096	58054.68146216384	3	0.0	740.4933859969877	0.0	0.0	14990.630786762184;	%6 OCGT
	2	989.8839972534882	395.9535989013953	3	0.0	193.7529843909744	0.0	0.0	4126.803481311164;	%6 biomass
	2	0.0	0.0	3	0.0	0.3251259869753945	0.0	0.0	5418.542086970701;	%6 onwind
	2	0.0	0.0	3	0.0	0.26544513975138007	0.0	0.0	1479.3554717745228;	%6 solar
	2	99606.87600339073	39842.750401356294	3	0.0	508.19834695607517	0.0	0.0	32918.80681525798;	%7 CCGT
	2	144516.21396244713	57806.48558497885	3	0.0	737.3276222573833	0.0	0.0	14990.630786762184;	%7 OCGT
	2	989.8829372216959	395.9531748886784	3	0.0	193.75277690775022	0.0	0.0	4126.803481311164;	%7 biomass
	2	0.0	0.0	3	0.0	0.32868219198720255	0.0	0.0	5418.542086970701;	%7 onwind
	2	0.0	0.0	3	0.0	0.25935135314048957	0.0	0.0	1479.355471774523;	%7 solar
	2	112619.65123488643	45047.86049395457	3	0.0	574.590057320849	0.0	0.0	32918.80681525798;	%8 CCGT
	2	147687.20553715067	59074.88221486027	3	0.0	753.5061506997483	0.0	0.0	14990.630786762184;	%8 OCGT
	2	989.8848189261455	395.9539275704582	3	0.0	193.75314521944517	0.0	0.0	4126.803481311164;	%8 biomass
	2	0.0	0.0	3	0.0	0.3322497716594224	0.0	0.0	5418.542086970701;	%8 onwind
	2	0.0	0.0	3	0.0	0.25902030549500354	0.0	0.0	1479.3554717745228;	%8 solar
	2	107992.69729309251	43197.078917237	3	0.0	550.9831494545535	0.0	0.0	32918.80681525798;	%9 CCGT
	2	145711.8833359165	58284.7533343666	3	0.0	743.4279762036556	0.0	0.0	14990.630786762184;	%9 OCGT
	2	989.8459705013673	395.93838820054697	3	0.0	193.7455412999349	0.0	0.0	4126.803481311164;	%9 biomass
	2	0.0	0.0	3	0.0	0.3231245183230015	0.0	0.0	5418.542086970701;	%9 onwind
	2	0.0	0.0	3	0.0	0.2620391260271471	0.0	0.0	1479.3554717745228;	%9 solar
	2	0.0	0.0	3	0.0	0.12568513004511467	0.0	0.0	263.1835406453565;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.094263847533836	0.0	0.0	263.1835406453565;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13693155306258936	0.0	0.0	263.1835406453565;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.10269866479694202	0.0	0.0	263.1835406453565;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1428807328835876	0.0	0.0	263.1835406453565;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.1071605496626907	0.0	0.0	263.1835406453565;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11931977324189913	0.0	0.0	263.1835406453565;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.08948982993142435	0.0	0.0	263.1835406453565;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1419559270884605	0.0	0.0	263.1835406453565;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.10646694531634537	0.0	0.0	263.1835406453565;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1365616979440249	0.0	0.0	263.1835406453565;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.10242127345801867	0.0	0.0	263.1835406453565;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12998860917488653	0.0	0.0	263.1835406453565;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.0974914568811649	0.0	0.0	263.1835406453565;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13194817976795656	0.0	0.0	263.1835406453565;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.13046305893629628	0.0	0.0	263.1835406453565;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.09784729420222221	0.0	0.0	263.1835406453565;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1270416592036898	0.0	0.0	263.1835406453565;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.09528124440276733	0.0	0.0	263.1835406453565;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12470659265828134	0.0	0.0	263.1835406453565;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.093529944493711	0.0	0.0	263.1835406453565;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.11942482877574509	0.0	0.0	263.1835406453565;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.08956862158180881	0.0	0.0	263.1835406453565;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13172187709441785	0.0	0.0	263.1835406453565;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.09879140782081339	0.0	0.0	263.1835406453565;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12091424700815223	0.0	0.0	263.1835406453565;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.09068568525611417	0.0	0.0	263.1835406453565;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12164169192308363	0.0	0.0	263.1835406453565;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.09123126894231273	0.0	0.0	263.1835406453565;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1409672561495313	0.0	0.0	263.1835406453565;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.10572544211214846	0.0	0.0	263.1835406453565;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12904458192978774	0.0	0.0	263.1835406453565;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.0967834364473408	0.0	0.0	263.1835406453565;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.13411700356563622	0.0	0.0	263.1835406453565;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.10058775267422716	0.0	0.0	263.1835406453565;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1359675495738037	0.0	0.0	263.1835406453565;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.10197566218035277	0.0	0.0	263.1835406453565;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.12072311301808634	0.0	0.0	263.1835406453565;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.1405204221031984	0.0	0.0	263.1835406453565;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.10539031657739881	0.0	0.0	263.1835406453565;	%DE1 76 PHS (pump mode)
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

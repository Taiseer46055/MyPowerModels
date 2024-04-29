function mpc = mpc_006_005
%MPC_006_005	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 6 	Weight: 24
%	Time step: 5

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 24;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1211.7692980620684	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	3779.207887276938	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	390.21680033131696	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-320.7473517972613	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	2413.9094165133433	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	4957.187008110202	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	185.9700019641561	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2152.826732587038	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2206.4796920653916	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	1603.313812933253	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1966.9126293065876	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	2893.557177024121	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4585.252825322475	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3909.068122541742	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	1426.4448226786972	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-2558.305164995927	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3036.429198881748	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3156.7850698511834	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5313.822472489543	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	2809.7347875534015	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8	8	469	10.0	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.580029381306893	0.0	0.0	1	1	276	17.0	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.726533190801673	0.0	0.0	43	43	124	38.0	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8491273597738445	0.0	0.0	23	23	168	28.0	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5	5	53	224.0	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2	2	53	224.0	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24	24	1172	10.0	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8151760257972884	0.0	0.0	57	57	309	38.0	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	75	75	419	28.0	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2	2	74	224.0	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1	1	74	224.0	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50	50	1637	10.0	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.286975836230926	0.0	0.0	20	20	963	17.0	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.331110456341534	0.0	0.0	118	118	431	38.0	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5341857499722489	0.0	0.0	86	86	585	28.0	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	46	224.0	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1	1	46	224.0	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35	35	1023	10.0	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.661103440015395	0.0	0.0	7	7	602	17.0	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	25.63451654741977	0.0	0.0	86	86	270	38.0	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.2084102666745442	0.0	0.0	40	40	366	28.0	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	58	224.0	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69	69	1280	10.0	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.25250951072135025	0.0	0.0	29	29	337	38.0	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.22777919368027072	0.0	0.0	166	166	457	28.0	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5	5	72	224.0	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3	3	72	224.0	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44	44	1597	10.0	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.786677730659278	0.0	0.0	45	45	421	38.0	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.11969889773487438	0.0	0.0	154	154	571	28.0	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1	1	28	224.0	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	28	224.0	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11	11	621	10.0	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.607698676625375	0.0	0.0	34	34	365	17.0	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.117172967977808	0.0	0.0	37	37	164	38.0	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.772973902911401	0.0	0.0	33	33	222	28.0	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2	2	75	224.0	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	75	224.0	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54	54	1672	10.0	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.229245910829034	0.0	0.0	46	46	984	17.0	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.473703337867185	0.0	0.0	137	137	440	38.0	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.2664333740563601	0.0	0.0	56	56	597	28.0	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2	2	58	224.0	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1	1	58	224.0	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52	52	1279	10.0	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.7207919339048401	0.0	0.0	18	18	337	38.0	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6919031723557258	0.0	0.0	173	173	457	28.0	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1	1	45	224.0	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1	1	45	224.0	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12	12	998	10.0	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.6918842440967996	0.0	0.0	72	72	263	38.0	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	67	67	357	28.0	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2	2	44	224.0	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	44	224.0	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18	18	974	10.0	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.635455933115752	0.0	0.0	53	53	257	38.0	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.013340847019512	0.0	0.0	81	81	348	28.0	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17	17	70	224.0	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4	4	70	224.0	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13	13	1562	10.0	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.643856832490269	0.0	0.0	49	49	411	38.0	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56	56	558	28.0	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4	4	81	224.0	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1	1	81	224.0	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63	63	1794	10.0	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.037331162746788	0.0	0.0	82	82	472	38.0	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.5709382848253407	0.0	0.0	155	155	641	28.0	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5	5	139	224.0	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1	1	139	224.0	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65	65	3109	10.0	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.380380447906068	0.0	0.0	350	350	1829	17.0	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.331559592031608	0.0	0.0	156	156	818	38.0	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.16243659623113008	0.0	0.0	73	73	1111	28.0	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12	12	58	224.0	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2	2	58	224.0	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38	38	1279	10.0	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.05437218967234381	0.0	0.0	16	16	337	38.0	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.17393458927170297	0.0	0.0	98	98	457	28.0	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	63	224.0	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1	1	63	224.0	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27	27	1391	10.0	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	17	17	366	38.0	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	90	90	497	28.0	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11	11	73	224.0	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1	1	73	224.0	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45	45	1622	10.0	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.488243496897965	0.0	0.0	70	70	427	38.0	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	91	91	580	28.0	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2	2	76	224.0	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2	2	76	224.0	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54	54	1682	10.0	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.294412994156561	0.0	0.0	109	109	443	38.0	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.2906548866211731	0.0	0.0	111	111	601	28.0	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5	5	101	224.0	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3	3	101	224.0	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60	60	2241	10.0	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.052006039256302	0.0	0.0	113	113	590	38.0	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.0007570611994236	0.0	0.0	170	170	801	28.0	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14	14	112	224.0	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1	1	112	224.0	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70	70	2491	10.0	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.265189114956547	0.0	0.0	133	133	656	38.0	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.12819957333438375	0.0	0.0	132	132	890	28.0	3;	%9 solar
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
	2	1827.4880191354994	730.9952076541997	3	0.0	357.6997492925229	0.0	0.0	7618.714119343689;	%0 biomass
	2	0.0	0.0	3	0.0	0.5921496559969553	0.0	0.0	8823.672871614346;	%0 offwind
	2	0.0	0.0	3	0.0	0.5783395734836745	0.0	0.0	10003.462314407449;	%0 onwind
	2	0.0	0.0	3	0.0	0.4574126569588455	0.0	0.0	2731.1177940452726;	%0 solar
	2	213189.66701372212	85275.86680548884	3	0.0	1087.7023827230719	0.0	0.0	60773.181812783965;	%1 CCGT
	2	260231.94410915227	104092.7776436609	3	0.0	1327.7140005568992	0.0	0.0	27675.010683253262;	%1 OCGT
	2	1827.4150146359184	730.9660058543675	3	0.0	357.68545990133464	0.0	0.0	7618.714119343687;	%1 biomass
	2	0.0	0.0	3	0.0	0.5899203358055567	0.0	0.0	10003.462314407449;	%1 onwind
	2	0.0	0.0	3	0.0	0.478554447350525	0.0	0.0	2731.1177940452726;	%1 solar
	2	186798.65172430128	74719.4606897205	3	0.0	953.0543455321492	0.0	0.0	60773.181812783965;	%10 CCGT
	2	270964.9785252792	108385.99141011169	3	0.0	1382.4743802310163	0.0	0.0	27675.010683253262;	%10 OCGT
	2	1827.4486317121464	730.9794526848584	3	0.0	357.6920398731936	0.0	0.0	7618.714119343687;	%10 biomass
	2	0.0	0.0	3	0.0	0.5926341249936754	0.0	0.0	8710.98731681927;	%10 offwind
	2	0.0	0.0	3	0.0	0.5951930224686903	0.0	0.0	10003.46231440745;	%10 onwind
	2	0.0	0.0	3	0.0	0.4757346722189284	0.0	0.0	2731.1177940452726;	%10 solar
	2	174301.92301128106	69720.76920451241	3	0.0	889.2955255677605	0.0	0.0	60773.181812783965;	%11 CCGT
	2	258751.31714702462	103500.52685880984	3	0.0	1320.1597813623703	0.0	0.0	27675.010683253262;	%11 OCGT
	2	1827.4890474401118	730.9956189760447	3	0.0	357.6999505656903	0.0	0.0	7618.714119343687;	%11 biomass
	2	0.0	0.0	3	0.0	0.6075413015623223	0.0	0.0	8661.12333920036;	%11 offwind
	2	0.0	0.0	3	0.0	0.5785634471459538	0.0	0.0	10003.462314407449;	%11 onwind
	2	0.0	0.0	3	0.0	0.48228966142461593	0.0	0.0	2731.1177940452726;	%11 solar
	2	257356.88715639265	102942.75486255705	3	0.0	1313.0453426346562	0.0	0.0	27675.010683253262;	%12 OCGT
	2	1827.4543148249554	730.9817259299822	3	0.0	357.6931522460277	0.0	0.0	7618.714119343687;	%12 biomass
	2	0.0	0.0	3	0.0	0.6003828699163838	0.0	0.0	10003.462314407449;	%12 onwind
	2	0.0	0.0	3	0.0	0.4802425738225203	0.0	0.0	2731.1177940452726;	%12 solar
	2	193484.4210130356	77393.76840521424	3	0.0	987.1654133318142	0.0	0.0	60773.181812783965;	%13 CCGT
	2	269145.3345992216	107658.13383968861	3	0.0	1373.1904826490895	0.0	0.0	27675.010683253262;	%13 OCGT
	2	1827.437539908096	730.9750159632383	3	0.0	357.6898688408878	0.0	0.0	7618.714119343687;	%13 biomass
	2	0.0	0.0	3	0.0	0.6103017298795039	0.0	0.0	10003.462314407449;	%13 onwind
	2	0.0	0.0	3	0.0	0.4858474120560355	0.0	0.0	2731.1177940452726;	%13 solar
	2	199384.70404655242	79753.88161862097	3	0.0	1017.268898196696	0.0	0.0	60773.181812783965;	%14 CCGT
	2	268952.6431494077	107581.05725976307	3	0.0	1372.207363007182	0.0	0.0	27675.010683253262;	%14 OCGT
	2	1827.362635861863	730.9450543447452	3	0.0	357.67520764569645	0.0	0.0	7618.714119343687;	%14 biomass
	2	0.0	0.0	3	0.0	0.5933930485307607	0.0	0.0	9121.698039774601;	%14 offwind
	2	0.0	0.0	3	0.0	0.5875200600520374	0.0	0.0	10003.46231440745;	%14 onwind
	2	0.0	0.0	3	0.0	0.4696767402684739	0.0	0.0	2731.1177940452726;	%14 solar
	2	179065.04409096672	71626.01763638669	3	0.0	913.5971637294219	0.0	0.0	60773.181812783965;	%15 CCGT
	2	251155.34446498242	100462.13778599298	3	0.0	1281.40481869889	0.0	0.0	27675.010683253262;	%15 OCGT
	2	1827.4583724529152	730.9833489811662	3	0.0	357.693946457803	0.0	0.0	7618.714119343687;	%15 biomass
	2	0.0	0.0	3	0.0	0.6007672661530072	0.0	0.0	9519.401770161521;	%15 offwind
	2	0.0	0.0	3	0.0	0.5946895270708468	0.0	0.0	10003.46231440745;	%15 onwind
	2	0.0	0.0	3	0.0	0.4641461664937011	0.0	0.0	2731.1177940452726;	%15 solar
	2	187440.78179122	74976.312716488	3	0.0	956.3305193429592	0.0	0.0	60773.181812783965;	%16 CCGT
	2	259347.27197201183	103738.90878880472	3	0.0	1323.200367204142	0.0	0.0	27675.010683253262;	%16 OCGT
	2	1827.4143366953808	730.9657346781523	3	0.0	357.68532720598563	0.0	0.0	7618.714119343689;	%16 biomass
	2	0.0	0.0	3	0.0	0.5936823166248113	0.0	0.0	10003.462314407449;	%16 onwind
	2	0.0	0.0	3	0.0	0.47741562038694174	0.0	0.0	2731.117794045273;	%16 solar
	2	192713.5574414266	77085.42297657064	3	0.0	983.2324359256457	0.0	0.0	60773.181812783965;	%17 CCGT
	2	252821.9728324652	101128.78913298607	3	0.0	1289.9080246554345	0.0	0.0	27675.010683253262;	%17 OCGT
	2	1827.4349732265007	730.9739892906002	3	0.0	357.6893664565474	0.0	0.0	7618.714119343687;	%17 biomass
	2	0.0	0.0	3	0.0	0.6118107353859376	0.0	0.0	10003.46231440745;	%17 onwind
	2	0.0	0.0	3	0.0	0.4892157951742765	0.0	0.0	2731.117794045273;	%17 solar
	2	202007.2411093918	80802.89644375672	3	0.0	1030.6491893336315	0.0	0.0	60773.181812783965;	%18 CCGT
	2	254789.58688153696	101915.83475261478	3	0.0	1299.9468718445762	0.0	0.0	27675.010683253262;	%18 OCGT
	2	1827.4680496280066	730.9872198512028	3	0.0	357.6958406005102	0.0	0.0	7618.714119343687;	%18 biomass
	2	0.0	0.0	3	0.0	0.5978302868628811	0.0	0.0	10003.462314407449;	%18 onwind
	2	0.0	0.0	3	0.0	0.47998642152732507	0.0	0.0	2731.1177940452726;	%18 solar
	2	189283.6028650283	75713.44114601133	3	0.0	965.7326676787159	0.0	0.0	60773.181812783965;	%19 CCGT
	2	265203.43513179483	106081.37405271793	3	0.0	1353.0787506724225	0.0	0.0	27675.010683253262;	%19 OCGT
	2	1827.4259412081133	730.9703764832454	3	0.0	357.68759859231034	0.0	0.0	7618.714119343687;	%19 biomass
	2	0.0	0.0	3	0.0	0.6099368606752811	0.0	0.0	10003.462314407449;	%19 onwind
	2	0.0	0.0	3	0.0	0.47709993180834764	0.0	0.0	2731.1177940452726;	%19 solar
	2	202645.96636440192	81058.38654576078	3	0.0	1033.9079916551118	0.0	0.0	60773.181812783965;	%2 CCGT
	2	267821.9153899025	107128.766155961	3	0.0	1366.438343826033	0.0	0.0	27675.010683253262;	%2 OCGT
	2	1827.41700402184	730.9668016087361	3	0.0	357.68584928984933	0.0	0.0	7618.714119343687;	%2 biomass
	2	0.0	0.0	3	0.0	0.6010887467605381	0.0	0.0	10003.462314407449;	%2 onwind
	2	0.0	0.0	3	0.0	0.4786324396376288	0.0	0.0	2731.1177940452726;	%2 solar
	2	220990.1173295057	88396.0469318023	3	0.0	1127.5005986199271	0.0	0.0	60773.181812783965;	%3 CCGT
	2	274895.2273383273	109958.09093533093	3	0.0	1402.5266700935067	0.0	0.0	27675.010683253262;	%3 OCGT
	2	1827.3636219316354	730.9454487726542	3	0.0	357.6754006521111	0.0	0.0	7618.714119343687;	%3 biomass
	2	0.0	0.0	3	0.0	0.5945833703578116	0.0	0.0	9477.018484767026;	%3 offwind
	2	0.0	0.0	3	0.0	0.611686009518218	0.0	0.0	10003.46231440745;	%3 onwind
	2	0.0	0.0	3	0.0	0.4792836574678448	0.0	0.0	2731.1177940452726;	%3 solar
	2	192639.96727885594	77055.98691154238	3	0.0	982.8569759125303	0.0	0.0	60773.181812783965;	%4 CCGT
	2	263451.59644486883	105380.63857794754	3	0.0	1344.1407981881064	0.0	0.0	27675.010683253262;	%4 OCGT
	2	1827.4721603082153	730.988864123286	3	0.0	357.69664519636234	0.0	0.0	7618.714119343687;	%4 biomass
	2	0.0	0.0	3	0.0	0.5983153830411345	0.0	0.0	10003.462314407449;	%4 onwind
	2	0.0	0.0	3	0.0	0.46906892771209696	0.0	0.0	2731.1177940452726;	%4 solar
	2	197866.35303026627	79146.5412121065	3	0.0	1009.5222093380931	0.0	0.0	60773.181812783965;	%5 CCGT
	2	258704.4780706957	103481.79122827828	3	0.0	1319.9208064831414	0.0	0.0	27675.010683253262;	%5 OCGT
	2	1827.4478460061187	730.9791384024475	3	0.0	357.6918860845799	0.0	0.0	7618.714119343687;	%5 biomass
	2	0.0	0.0	3	0.0	0.5910414570712378	0.0	0.0	10003.462314407449;	%5 onwind
	2	0.0	0.0	3	0.0	0.48040813099006746	0.0	0.0	2731.1177940452726;	%5 solar
	2	194006.0169883385	77602.40679533541	3	0.0	989.8266172874414	0.0	0.0	60773.181812783965;	%6 CCGT
	2	267944.68367152545	107177.87346861017	3	0.0	1367.0647126098236	0.0	0.0	27675.010683253262;	%6 OCGT
	2	1827.4781487756707	730.9912595102683	3	0.0	357.69781733718355	0.0	0.0	7618.714119343687;	%6 biomass
	2	0.0	0.0	3	0.0	0.6002325913391898	0.0	0.0	10003.462314407449;	%6 onwind
	2	0.0	0.0	3	0.0	0.49005256569485545	0.0	0.0	2731.1177940452726;	%6 solar
	2	183889.61723702904	73555.84689481162	3	0.0	938.2123328419849	0.0	0.0	60773.181812783965;	%7 CCGT
	2	266799.16423836397	106719.66569534557	3	0.0	1361.2202257059382	0.0	0.0	27675.010683253262;	%7 OCGT
	2	1827.4761917939002	730.9904767175601	3	0.0	357.69743429123116	0.0	0.0	7618.714119343687;	%7 biomass
	2	0.0	0.0	3	0.0	0.6067978928994509	0.0	0.0	10003.462314407449;	%7 onwind
	2	0.0	0.0	3	0.0	0.47880249810551917	0.0	0.0	2731.117794045273;	%7 solar
	2	207913.20227979033	83165.28091191612	3	0.0	1060.7816442846442	0.0	0.0	60773.181812783965;	%8 CCGT
	2	272653.3025301243	109061.32101204974	3	0.0	1391.0882782149201	0.0	0.0	27675.010683253262;	%8 OCGT
	2	1827.479665709807	730.9918662839227	3	0.0	357.6981142512834	0.0	0.0	7618.714119343687;	%8 biomass
	2	0.0	0.0	3	0.0	0.6133841938327798	0.0	0.0	10003.462314407449;	%8 onwind
	2	0.0	0.0	3	0.0	0.47819133322154495	0.0	0.0	2731.1177940452726;	%8 solar
	2	199371.1334641708	79748.45338566831	3	0.0	1017.1996605314835	0.0	0.0	60773.181812783965;	%9 CCGT
	2	269006.55385092273	107602.62154036912	3	0.0	1372.4824176067486	0.0	0.0	27675.010683253262;	%9 OCGT
	2	1827.4079455409858	730.9631782163943	3	0.0	357.68407624603367	0.0	0.0	7618.714119343687;	%9 biomass
	2	0.0	0.0	3	0.0	0.5965375722886181	0.0	0.0	10003.462314407449;	%9 onwind
	2	0.0	0.0	3	0.0	0.48376454035780997	0.0	0.0	2731.1177940452726;	%9 solar
	2	0.0	0.0	3	0.0	0.23203408623713478	0.0	0.0	485.87730580681205;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.17402556467785107	0.0	0.0	485.87730580681205;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2527967133463188	0.0	0.0	485.87730580681205;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.1895975350097391	0.0	0.0	485.87730580681205;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26377981455431554	0.0	0.0	485.87730580681205;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.19783486091573665	0.0	0.0	485.87730580681205;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22028265829273685	0.0	0.0	485.87730580681205;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.16521199371955264	0.0	0.0	485.87730580681205;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.26207248077869627	0.0	0.0	485.87730580681205;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.1965543605840222	0.0	0.0	485.87730580681205;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25211390389666133	0.0	0.0	485.87730580681205;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.189085427922496	0.0	0.0	485.87730580681205;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23997897078440592	0.0	0.0	485.87730580681205;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.17998422808830444	0.0	0.0	485.87730580681205;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2435966395716121	0.0	0.0	485.87730580681205;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.2408548780362393	0.0	0.0	485.87730580681205;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.1806411585271795	0.0	0.0	485.87730580681205;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23453844776065808	0.0	0.0	485.87730580681205;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.17590383582049357	0.0	0.0	485.87730580681205;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23022755567682707	0.0	0.0	485.87730580681205;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.1726706667576203	0.0	0.0	485.87730580681205;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22047660697060634	0.0	0.0	485.87730580681205;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.16535745522795475	0.0	0.0	485.87730580681205;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.24317885002046372	0.0	0.0	485.87730580681205;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.1823841375153478	0.0	0.0	485.87730580681205;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22322630216889644	0.0	0.0	485.87730580681205;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.16741972662667232	0.0	0.0	485.87730580681205;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2245692773964621	0.0	0.0	485.87730580681205;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.1684269580473466	0.0	0.0	485.87730580681205;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2602472421222116	0.0	0.0	485.87730580681205;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.19518543159165871	0.0	0.0	485.87730580681205;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.23823615125499276	0.0	0.0	485.87730580681205;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.17867711344124457	0.0	0.0	485.87730580681205;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2476006219673284	0.0	0.0	485.87730580681205;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.18570046647549632	0.0	0.0	485.87730580681205;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.25101701459779147	0.0	0.0	485.87730580681205;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.1882627609483436	0.0	0.0	485.87730580681205;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.22287343941800555	0.0	0.0	485.87730580681205;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.25942231772898167	0.0	0.0	485.87730580681205;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.19456673829673626	0.0	0.0	485.87730580681205;	%DE1 76 PHS (pump mode)
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

function mpc = mpc_001_012
%MPC_001_012	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 1 	Weight: 40
%	Time step: 12

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 40;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	807.7070356094348	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5863.754933940212	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1278.6246216523098	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-240.41832828361785	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	1145.1054462464044	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	7116.754220311652	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	294.0710842455114	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2977.796769634356	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2252.9423413143513	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	4748.937746055479	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	215.3487774911082	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	5411.039025144165	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5864.005022781655	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	4902.042026713038	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3612.67564979023	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	2963.9048550806356	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	7914.26380818223	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5307.397078319469	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4468.158187866036	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4320.156827607622	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8.0	8.0	469.0	10	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.791387538174586	0.0	0.0	1.0	1.0	276.0	17	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.675536974198153	0.0	0.0	43.0	43.0	124.0	38	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.328087168928151	0.0	0.0	23.0	23.0	168.0	28	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5.0	5.0	53.0	224	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2.0	2.0	53.0	224	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24.0	24.0	1172.0	10	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.335480129051924	0.0	0.0	57.0	57.0	309.0	38	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.083752378549571	0.0	0.0	75.0	75.0	419.0	28	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2.0	2.0	74.0	224	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1.0	1.0	74.0	224	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50.0	50.0	1637.0	10	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.43551413959386	0.0	0.0	20.0	20.0	963.0	17	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.287883854740299	0.0	0.0	118.0	118.0	431.0	38	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.525700836230206	0.0	0.0	86.0	86.0	585.0	28	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	46.0	224	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1.0	1.0	46.0	224	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35.0	35.0	1023.0	10	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.065253879792355	0.0	0.0	7.0	7.0	602.0	17	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	20.600161365327963	0.0	0.0	86.0	86.0	270.0	38	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.775188061041026	0.0	0.0	40.0	40.0	366.0	28	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	58.0	224	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69.0	69.0	1280.0	10	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0855101704738113	0.0	0.0	29.0	29.0	337.0	38	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.377326101755347	0.0	0.0	166.0	166.0	457.0	28	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5.0	5.0	72.0	224	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3.0	3.0	72.0	224	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44.0	44.0	1597.0	10	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1562097556751416	0.0	0.0	45.0	45.0	421.0	38	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.04597762059597	0.0	0.0	154.0	154.0	571.0	28	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	28.0	224	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	28.0	224	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11.0	11.0	621.0	10	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.63364415392492	0.0	0.0	34.0	34.0	365.0	17	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	21.155846967871966	0.0	0.0	37.0	37.0	164.0	38	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.303105350732423	0.0	0.0	33.0	33.0	222.0	28	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2.0	2.0	75.0	224	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	75.0	224	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54.0	54.0	1672.0	10	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.785340361074777	0.0	0.0	46.0	46.0	984.0	17	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.636002554971522	0.0	0.0	137.0	137.0	440.0	38	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.988378770548093	0.0	0.0	56.0	56.0	597.0	28	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2.0	2.0	58.0	224	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1.0	1.0	58.0	224	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52.0	52.0	1279.0	10	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.344917590853193	0.0	0.0	18.0	18.0	337.0	38	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.739681187578994	0.0	0.0	173.0	173.0	457.0	28	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1.0	1.0	45.0	224	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1.0	1.0	45.0	224	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12.0	12.0	991.0	10	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.299133580700223	0.0	0.0	72.0	72.0	261.0	38	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.342607546704661	0.0	0.0	67.0	67.0	354.0	28	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2.0	2.0	44.0	224	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	44.0	224	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18.0	18.0	974.0	10	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.77663411877554	0.0	0.0	53.0	53.0	257.0	38	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.831589699264869	0.0	0.0	81.0	81.0	348.0	28	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17.0	17.0	70.0	224	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4.0	4.0	70.0	224	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13.0	13.0	1562.0	10	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.737236422203777	0.0	0.0	49.0	49.0	411.0	38	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.488638342208655	0.0	0.0	56.0	56.0	558.0	28	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4.0	4.0	81.0	224	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1.0	1.0	81.0	224	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63.0	63.0	1794.0	10	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.652061356398702	0.0	0.0	82.0	82.0	472.0	38	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.374983871643713	0.0	0.0	155.0	155.0	641.0	28	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5.0	5.0	139.0	224	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	139.0	224	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65.0	65.0	3109.0	10	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.07914734153271	0.0	0.0	350.0	350.0	1829.0	17	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.55564277769378	0.0	0.0	156.0	156.0	818.0	38	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.714003888468828	0.0	0.0	73.0	73.0	1111.0	28	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12.0	12.0	58.0	224	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2.0	2.0	58.0	224	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38.0	38.0	1279.0	10	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.329360149225998	0.0	0.0	16.0	16.0	337.0	38	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.06243772592433	0.0	0.0	98.0	98.0	457.0	28	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2.0	2.0	35.0	224	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1.0	1.0	35.0	224	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27.0	27.0	779.0	10	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.26217020725442	0.0	0.0	17.0	17.0	205.0	38	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.746913477426173	0.0	0.0	90.0	90.0	278.0	28	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11.0	11.0	73.0	224	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1.0	1.0	73.0	224	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45.0	45.0	1622.0	10	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.609690671597993	0.0	0.0	70.0	70.0	427.0	38	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.759155848046736	0.0	0.0	91.0	91.0	580.0	28	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2.0	2.0	76.0	224	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26000000000002	78.4	2.8	2.0	2.0	76.0	224	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54.0	54.0	1682.0	10	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.2156014066928025	0.0	0.0	109.0	109.0	443.0	38	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.254061487939173	0.0	0.0	111.0	111.0	601.0	28	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5.0	5.0	101.0	224	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3.0	3.0	101.0	224	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60.0	60.0	2241.0	10	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.339919697139504	0.0	0.0	113.0	113.0	590.0	38	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.5647929346665155	0.0	0.0	170.0	170.0	801.0	28	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14.0	14.0	112.0	224	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1.0	1.0	112.0	224	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70.0	70.0	2491.0	10	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.173651943278252	0.0	0.0	133.0	133.0	656.0	38	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.274379315722102	0.0	0.0	132.0	132.0	890.0	28	3;	%9 solar
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
	2	3045.813365225832	1218.325346090333	3	0.0	596.1662488208715	0.0	0.0	11123322.614241786;	%0 biomass
	2	0.0	0.0	3	0.0	0.9869160933282589	0.0	0.0	7577977.877974673;	%0 offwind
	2	0.0	0.0	3	0.0	0.9638992891394573	0.0	0.0	3843435.5207986515;	%0 onwind
	2	0.0	0.0	3	0.0	0.7623544282647425	0.0	0.0	1424082.849752178;	%0 solar
	2	355316.11168953683	142126.44467581474	3	0.0	1812.8373045384528	0.0	0.0	3961109.171726098;	%1 CCGT
	2	433719.9068485871	173487.96273943485	3	0.0	2212.8566675948323	0.0	0.0	1803817.6606048997;	%1 OCGT
	2	3045.6916910598643	1218.2766764239457	3	0.0	596.142433168891	0.0	0.0	11123322.614241784;	%1 biomass
	2	0.0	0.0	3	0.0	0.9832005596759277	0.0	0.0	3843435.5207986515;	%1 onwind
	2	0.0	0.0	3	0.0	0.7975907455842084	0.0	0.0	1424082.849752178;	%1 solar
	2	311331.08620716876	124532.43448286751	3	0.0	1588.4239092202486	0.0	0.0	3961109.171726098;	%10 CCGT
	2	451608.297542132	180643.3190168528	3	0.0	2304.123967051694	0.0	0.0	1803817.6606048997;	%10 OCGT
	2	3045.7477195202437	1218.2990878080975	3	0.0	596.153399788656	0.0	0.0	11123322.614241784;	%10 biomass
	2	0.0	0.0	3	0.0	0.9877235416561257	0.0	0.0	7481200.872091845;	%10 offwind
	2	0.0	0.0	3	0.0	0.9919883707811505	0.0	0.0	3843435.5207986524;	%10 onwind
	2	0.0	0.0	3	0.0	0.7928911203648807	0.0	0.0	1424082.849752178;	%10 solar
	2	290503.20501880173	116201.2820075207	3	0.0	1482.1592092796006	0.0	0.0	3961109.171726098;	%11 CCGT
	2	431252.195245041	172500.8780980164	3	0.0	2200.266302270617	0.0	0.0	1803817.6606048997;	%11 OCGT
	2	3045.8150790668533	1218.3260316267413	3	0.0	596.1665842761505	0.0	0.0	11123322.614241783;	%11 biomass
	2	0.0	0.0	3	0.0	1.012568835937204	0.0	0.0	7438376.514842661;	%11 offwind
	2	0.0	0.0	3	0.0	0.9642724119099229	0.0	0.0	3843435.520798652;	%11 onwind
	2	0.0	0.0	3	0.0	0.80381610237436	0.0	0.0	1424082.849752178;	%11 solar
	2	428928.1452606544	171571.25810426177	3	0.0	2188.4089043910935	0.0	0.0	1803817.6606048997;	%12 OCGT
	2	3045.7571913749257	1218.3028765499703	3	0.0	596.1552537433795	0.0	0.0	11123322.614241784;	%12 biomass
	2	0.0	0.0	3	0.0	1.0006381165273062	0.0	0.0	3843435.520798652;	%12 onwind
	2	0.0	0.0	3	0.0	0.8004042897042005	0.0	0.0	1424082.849752178;	%12 solar
	2	322474.03502172604	128989.6140086904	3	0.0	1645.2756888863569	0.0	0.0	3961109.171726098;	%13 CCGT
	2	448575.55766536924	179430.2230661477	3	0.0	2288.650804415149	0.0	0.0	1803817.6606048997;	%13 OCGT
	2	3045.7292331801596	1218.2916932720639	3	0.0	596.1497814014797	0.0	0.0	11123322.614241784;	%13 biomass
	2	0.0	0.0	3	0.0	1.017169549799173	0.0	0.0	3843435.520798652;	%13 onwind
	2	0.0	0.0	3	0.0	0.8097456867600592	0.0	0.0	1424082.849752178;	%13 solar
	2	332307.8400775874	132923.13603103496	3	0.0	1695.4481636611602	0.0	0.0	3961109.171726098;	%14 CCGT
	2	448254.40524901275	179301.7620996051	3	0.0	2287.012271678637	0.0	0.0	1803817.6606048997;	%14 OCGT
	2	3045.6043931031054	1218.241757241242	3	0.0	596.1253460761607	0.0	0.0	11123322.614241784;	%14 biomass
	2	0.0	0.0	3	0.0	0.9889884142179345	0.0	0.0	7833928.9047476;	%14 offwind
	2	0.0	0.0	3	0.0	0.979200100086729	0.0	0.0	3843435.5207986524;	%14 onwind
	2	0.0	0.0	3	0.0	0.7827945671141232	0.0	0.0	1424082.849752178;	%14 solar
	2	298441.7401516112	119376.69606064448	3	0.0	1522.6619395490366	0.0	0.0	3961109.171726098;	%15 CCGT
	2	418592.24077497073	167436.89630998828	3	0.0	2135.674697831483	0.0	0.0	1803817.6606048997;	%15 OCGT
	2	3045.7639540881923	1218.305581635277	3	0.0	596.1565774296716	0.0	0.0	11123322.614241784;	%15 biomass
	2	0.0	0.0	3	0.0	1.0012787769216787	0.0	0.0	8175486.226138718;	%15 offwind
	2	0.0	0.0	3	0.0	0.9911492117847447	0.0	0.0	3843435.5207986524;	%15 onwind
	2	0.0	0.0	3	0.0	0.7735769441561685	0.0	0.0	1424082.849752178;	%15 solar
	2	312401.3029853667	124960.52119414668	3	0.0	1593.884198904932	0.0	0.0	3961109.171726098;	%16 CCGT
	2	432245.4532866864	172898.18131467456	3	0.0	2205.3339453402364	0.0	0.0	1803817.6606049002;	%16 OCGT
	2	3045.6905611589677	1218.2762244635871	3	0.0	596.142212009976	0.0	0.0	11123322.614241786;	%16 biomass
	2	0.0	0.0	3	0.0	0.9894705277080189	0.0	0.0	3843435.520798652;	%16 onwind
	2	0.0	0.0	3	0.0	0.7956927006449028	0.0	0.0	1424082.8497521784;	%16 solar
	2	321189.2624023777	128475.70496095106	3	0.0	1638.7207265427428	0.0	0.0	3961109.171726098;	%17 CCGT
	2	421369.95472077536	168547.98188831014	3	0.0	2149.8467077590576	0.0	0.0	1803817.6606048997;	%17 OCGT
	2	3045.724955377501	1218.2899821510005	3	0.0	596.1489440942456	0.0	0.0	11123322.614241784;	%17 biomass
	2	0.0	0.0	3	0.0	1.0196845589765626	0.0	0.0	3843435.5207986524;	%17 onwind
	2	0.0	0.0	3	0.0	0.8153596586237942	0.0	0.0	1424082.8497521784;	%17 solar
	2	336678.73518231965	134671.49407292786	3	0.0	1717.748648889386	0.0	0.0	3961109.171726098;	%18 CCGT
	2	424649.31146922824	169859.7245876913	3	0.0	2166.5781197409606	0.0	0.0	1803817.6606048997;	%18 OCGT
	2	3045.7800827133447	1218.3120330853378	3	0.0	596.1597343341837	0.0	0.0	11123322.614241784;	%18 biomass
	2	0.0	0.0	3	0.0	0.9963838114381351	0.0	0.0	3843435.520798652;	%18 onwind
	2	0.0	0.0	3	0.0	0.7999773692122085	0.0	0.0	1424082.849752178;	%18 solar
	2	315472.6714417139	126189.06857668555	3	0.0	1609.5544461311931	0.0	0.0	3961109.171726098;	%19 CCGT
	2	442005.72521965805	176802.29008786322	3	0.0	2255.1312511207043	0.0	0.0	1803817.6606048997;	%19 OCGT
	2	3045.7099020135224	1218.283960805409	3	0.0	596.1459976538505	0.0	0.0	11123322.614241784;	%19 biomass
	2	0.0	0.0	3	0.0	1.0165614344588016	0.0	0.0	3843435.5207986515;	%19 onwind
	2	0.0	0.0	3	0.0	0.7951665530139128	0.0	0.0	1424082.849752178;	%19 solar
	2	337743.2772740032	135097.31090960128	3	0.0	1723.1799860918532	0.0	0.0	3961109.171726098;	%2 CCGT
	2	446369.85898317076	178547.9435932683	3	0.0	2277.3972397100547	0.0	0.0	1803817.6606049002;	%2 OCGT
	2	3045.695006703067	1218.2780026812268	3	0.0	596.1430821497489	0.0	0.0	11123322.614241784;	%2 biomass
	2	0.0	0.0	3	0.0	1.0018145779342302	0.0	0.0	3843435.520798652;	%2 onwind
	2	0.0	0.0	3	0.0	0.7977207327293813	0.0	0.0	1424082.849752178;	%2 solar
	2	368316.8622158428	147326.74488633714	3	0.0	1879.167664366545	0.0	0.0	3961109.171726098;	%3 CCGT
	2	458158.7122305456	183263.48489221823	3	0.0	2337.5444501558445	0.0	0.0	1803817.6606048997;	%3 OCGT
	2	3045.606036552726	1218.2424146210903	3	0.0	596.1256677535185	0.0	0.0	11123322.614241784;	%3 biomass
	2	0.0	0.0	3	0.0	0.9909722839296862	0.0	0.0	8139086.463388152;	%3 offwind
	2	0.0	0.0	3	0.0	1.0194766825303634	0.0	0.0	3843435.5207986524;	%3 onwind
	2	0.0	0.0	3	0.0	0.7988060957797414	0.0	0.0	1424082.849752178;	%3 solar
	2	321066.6121314266	128426.64485257062	3	0.0	1638.094959854217	0.0	0.0	3961109.171726098;	%4 CCGT
	2	439085.9940747814	175634.39762991256	3	0.0	2240.2346636468437	0.0	0.0	1803817.6606048997;	%4 OCGT
	2	3045.7869338470255	1218.3147735388102	3	0.0	596.1610753272705	0.0	0.0	11123322.614241784;	%4 biomass
	2	0.0	0.0	3	0.0	0.9971923050685576	0.0	0.0	3843435.520798652;	%4 onwind
	2	0.0	0.0	3	0.0	0.7817815461868283	0.0	0.0	1424082.849752178;	%4 solar
	2	329777.2550504438	131910.9020201775	3	0.0	1682.5370155634887	0.0	0.0	3961109.171726098;	%5 CCGT
	2	431174.1301178262	172469.65204713048	3	0.0	2199.8680108052354	0.0	0.0	1803817.6606048997;	%5 OCGT
	2	3045.7464100101975	1218.298564004079	3	0.0	596.1531434742998	0.0	0.0	11123322.614241783;	%5 biomass
	2	0.0	0.0	3	0.0	0.9850690951187296	0.0	0.0	3843435.520798652;	%5 onwind
	2	0.0	0.0	3	0.0	0.8006802183167792	0.0	0.0	1424082.849752178;	%5 solar
	2	323343.36164723086	129337.34465889234	3	0.0	1649.7110288124022	0.0	0.0	3961109.171726098;	%6 CCGT
	2	446574.47278587567	178629.78911435028	3	0.0	2278.441187683039	0.0	0.0	1803817.6606049002;	%6 OCGT
	2	3045.7969146261175	1218.318765850447	3	0.0	596.1630288953058	0.0	0.0	11123322.614241784;	%6 biomass
	2	0.0	0.0	3	0.0	1.000387652231983	0.0	0.0	3843435.520798652;	%6 onwind
	2	0.0	0.0	3	0.0	0.8167542761580925	0.0	0.0	1424082.849752178;	%6 solar
	2	306482.6953950484	122593.07815801936	3	0.0	1563.6872214033083	0.0	0.0	3961109.171726098;	%7 CCGT
	2	444665.27373060654	177866.10949224263	3	0.0	2268.7003761765636	0.0	0.0	1803817.6606048997;	%7 OCGT
	2	3045.7936529898334	1218.3174611959334	3	0.0	596.1623904853853	0.0	0.0	11123322.614241783;	%7 biomass
	2	0.0	0.0	3	0.0	1.0113298214990847	0.0	0.0	3843435.520798652;	%7 onwind
	2	0.0	0.0	3	0.0	0.7980041635091986	0.0	0.0	1424082.8497521784;	%7 solar
	2	346522.0037996506	138608.80151986022	3	0.0	1767.969407141074	0.0	0.0	3961109.171726098;	%8 CCGT
	2	454422.17088354053	181768.8683534162	3	0.0	2318.4804636915333	0.0	0.0	1803817.6606048997;	%8 OCGT
	2	3045.7994428496786	1218.3197771398713	3	0.0	596.163523752139	0.0	0.0	11123322.614241784;	%8 biomass
	2	0.0	0.0	3	0.0	1.0223069897212995	0.0	0.0	3843435.520798652;	%8 onwind
	2	0.0	0.0	3	0.0	0.7969855553692416	0.0	0.0	1424082.849752178;	%8 solar
	2	332285.2224402846	132914.08897611385	3	0.0	1695.3327675524724	0.0	0.0	3961109.171726098;	%9 CCGT
	2	448344.25641820463	179337.70256728184	3	0.0	2287.470696011248	0.0	0.0	1803817.6606048997;	%9 OCGT
	2	3045.6799092349765	1218.2719636939905	3	0.0	596.1401270767227	0.0	0.0	11123322.614241784;	%9 biomass
	2	0.0	0.0	3	0.0	0.994229287147697	0.0	0.0	3843435.5207986515;	%9 onwind
	2	0.0	0.0	3	0.0	0.8062742339296833	0.0	0.0	1424082.849752178;	%9 solar
	2	0.0	0.0	3	0.0	0.3867234770618913	0.0	0.0	7093808.664779456;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.2900426077964185	0.0	0.0	7093808.664779456;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.421327855577198	0.0	0.0	7093808.664779456;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.3159958916828985	0.0	0.0	7093808.664779456;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43963302425719264	0.0	0.0	7093808.664779456;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.3297247681928945	0.0	0.0	7093808.664779456;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3671377638212281	0.0	0.0	7093808.664779456;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.27535332286592107	0.0	0.0	7093808.664779456;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43678746796449386	0.0	0.0	7093808.664779456;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.3275906009733704	0.0	0.0	7093808.664779456;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4201898398277689	0.0	0.0	7093808.664779456;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.3151423798708267	0.0	0.0	7093808.664779456;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3999649513073432	0.0	0.0	7093808.664779456;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.2999737134805074	0.0	0.0	7093808.664779456;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40599439928602016	0.0	0.0	7093808.664779456;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.40142479672706555	0.0	0.0	7093808.664779456;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.30106859754529913	0.0	0.0	7093808.664779456;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3908974129344301	0.0	0.0	7093808.664779456;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.2931730597008226	0.0	0.0	7093808.664779456;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3837125927947118	0.0	0.0	7093808.664779456;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.28778444459603386	0.0	0.0	7093808.664779456;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3674610116176772	0.0	0.0	7093808.664779456;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.2755957587132579	0.0	0.0	7093808.664779456;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40529808336743955	0.0	0.0	7093808.664779456;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.30397356252557967	0.0	0.0	7093808.664779456;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37204383694816073	0.0	0.0	7093808.664779456;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.27903287771112056	0.0	0.0	7093808.664779456;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37428212899410346	0.0	0.0	7093808.664779456;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.2807115967455776	0.0	0.0	7093808.664779456;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4337454035370193	0.0	0.0	7093808.664779456;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.32530905265276444	0.0	0.0	7093808.664779456;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39706025209165463	0.0	0.0	7093808.664779456;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.297795189068741	0.0	0.0	7093808.664779456;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41266770327888064	0.0	0.0	7093808.664779456;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.3095007774591605	0.0	0.0	7093808.664779456;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41836169099631904	0.0	0.0	7093808.664779456;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.3137712682472393	0.0	0.0	7093808.664779456;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3714557323633426	0.0	0.0	7093808.664779456;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.43237052954830274	0.0	0.0	7093808.664779456;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.32427789716122707	0.0	0.0	7093808.664779456;	%DE1 76 PHS (pump mode)
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

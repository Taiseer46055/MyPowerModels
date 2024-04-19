function mpc = mpc_003_007
%MPC_003_007	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 3 	Weight: 60
%	Time step: 7

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 60;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	366.98113685567324	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4449.320285556475	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1312.4904786199647	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	713.6049878655776	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	-432.4089300582554	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	5840.334285375703	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	220.63166460698818	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2506.9398781810455	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	2424.2048579043403	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	3720.793011388686	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1009.5919866472143	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	5349.4055757175065	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5166.282102888753	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3731.8860612737203	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	1405.8195008708053	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	197.10755426284825	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	7075.840846242953	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4263.4779487761925	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3670.9457102441693	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3395.7341340397275	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8.0	8.0	469.0	10	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.0785294548937983	0.0	0.0	1.0	1.0	276.0	17	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.876508354307768	0.0	0.0	43.0	43.0	124.0	38	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.507776218750184	0.0	0.0	23.0	23.0	168.0	28	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5.0	5.0	53.0	224	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2.0	2.0	53.0	224	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24.0	24.0	1172.0	10	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8848402474845536	0.0	0.0	57.0	57.0	309.0	38	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.536966097349508	0.0	0.0	75.0	75.0	419.0	28	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2.0	2.0	74.0	224	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1.0	1.0	74.0	224	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50.0	50.0	1637.0	10	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.130774429884465	0.0	0.0	20.0	20.0	963.0	17	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.42125866435775	0.0	0.0	118.0	118.0	431.0	38	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.087894926277494	0.0	0.0	86.0	86.0	585.0	28	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	46.0	224	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1.0	1.0	46.0	224	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35.0	35.0	1023.0	10	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.035224763903733	0.0	0.0	7.0	7.0	602.0	17	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.345728925341127	0.0	0.0	86.0	86.0	270.0	38	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.823050478582563	0.0	0.0	40.0	40.0	366.0	28	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	58.0	224	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69.0	69.0	1280.0	10	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.7681103482441278	0.0	0.0	29.0	29.0	337.0	38	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.061806572401892	0.0	0.0	166.0	166.0	457.0	28	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.29200000000003	78.4	5.0	5.0	5.0	72.0	224	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3.0	3.0	72.0	224	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44.0	44.0	1597.0	10	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1291192446426392	0.0	0.0	45.0	45.0	421.0	38	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.814560884591668	0.0	0.0	154.0	154.0	571.0	28	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	28.0	224	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	28.0	224	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11.0	11.0	621.0	10	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6573702325318544	0.0	0.0	34.0	34.0	365.0	17	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.672153273054359	0.0	0.0	37.0	37.0	164.0	38	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.110720071817655	0.0	0.0	33.0	33.0	222.0	28	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2.0	2.0	75.0	224	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	75.0	224	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54.0	54.0	1672.0	10	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3698300015746003	0.0	0.0	46.0	46.0	984.0	17	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.512012988710602	0.0	0.0	137.0	137.0	440.0	38	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.806570728078341	0.0	0.0	56.0	56.0	597.0	28	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2.0	2.0	58.0	224	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1.0	1.0	58.0	224	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52.0	52.0	1279.0	10	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.5885325102974914	0.0	0.0	18.0	18.0	337.0	38	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.731715529046195	0.0	0.0	173.0	173.0	457.0	28	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1.0	1.0	45.0	224	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1.0	1.0	45.0	224	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12.0	12.0	991.0	10	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.5314430418705727	0.0	0.0	72.0	72.0	261.0	38	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.912915761324751	0.0	0.0	67.0	67.0	354.0	28	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2.0	2.0	44.0	224	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	44.0	224	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18.0	18.0	974.0	10	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.092323342691421	0.0	0.0	53.0	53.0	257.0	38	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.656254439228086	0.0	0.0	81.0	81.0	348.0	28	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17.0	17.0	70.0	224	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4.0	4.0	70.0	224	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13.0	13.0	1562.0	10	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6040018842595796	0.0	0.0	49.0	49.0	411.0	38	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.713693838165442	0.0	0.0	56.0	56.0	558.0	28	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4.0	4.0	81.0	224	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1.0	1.0	81.0	224	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63.0	63.0	1794.0	10	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1981558553145693	0.0	0.0	82.0	82.0	472.0	38	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.180317825961847	0.0	0.0	155.0	155.0	641.0	28	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5.0	5.0	139.0	224	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	139.0	224	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65.0	65.0	3109.0	10	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.162247220251439	0.0	0.0	350.0	350.0	1829.0	17	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.614008084347549	0.0	0.0	156.0	156.0	818.0	38	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.567199309101824	0.0	0.0	73.0	73.0	1111.0	28	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12.0	12.0	58.0	224	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2.0	2.0	58.0	224	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38.0	38.0	1279.0	10	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.180664456932467	0.0	0.0	16.0	16.0	337.0	38	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.499497658047408	0.0	0.0	98.0	98.0	457.0	28	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2.0	2.0	35.0	224	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1.0	1.0	35.0	224	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27.0	27.0	779.0	10	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8069246414756421	0.0	0.0	17.0	17.0	205.0	38	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.210006610660024	0.0	0.0	90.0	90.0	278.0	28	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11.0	11.0	73.0	224	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1.0	1.0	73.0	224	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45.0	45.0	1622.0	10	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.5386614233022557	0.0	0.0	70.0	70.0	427.0	38	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.71710011806649	0.0	0.0	91.0	91.0	580.0	28	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2.0	2.0	76.0	224	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2.0	2.0	76.0	224	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54.0	54.0	1682.0	10	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1377616786997473	0.0	0.0	109.0	109.0	443.0	38	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.854595413653573	0.0	0.0	111.0	111.0	601.0	28	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5.0	5.0	101.0	224	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3.0	3.0	101.0	224	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60.0	60.0	2241.0	10	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.060599402147104	0.0	0.0	113.0	113.0	590.0	38	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.32283649104482	0.0	0.0	170.0	170.0	801.0	28	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14.0	14.0	112.0	224	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1.0	1.0	112.0	224	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70.0	70.0	2491.0	10	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.29678583737959	0.0	0.0	133.0	133.0	656.0	38	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.079771365777885	0.0	0.0	132.0	132.0	890.0	28	3;	%9 solar
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
	2	4568.7200478387485	1827.4880191354994	3	0.0	894.2493732313071	0.0	0.0	16684983.921362681;	%0 biomass
	2	0.0	0.0	3	0.0	1.4803741399923884	0.0	0.0	11366966.81696201;	%0 offwind
	2	0.0	0.0	3	0.0	1.445848933709186	0.0	0.0	5765153.281197977;	%0 onwind
	2	0.0	0.0	3	0.0	1.1435316423971138	0.0	0.0	2136124.274628267;	%0 solar
	2	532974.1675343052	213189.6670137221	3	0.0	2719.2559568076795	0.0	0.0	5941663.7575891465;	%1 CCGT
	2	650579.8602728806	260231.94410915227	3	0.0	3319.285001392248	0.0	0.0	2705726.4909073496;	%1 OCGT
	2	4568.537536589796	1827.4150146359186	3	0.0	894.2136497533365	0.0	0.0	16684983.921362678;	%1 biomass
	2	0.0	0.0	3	0.0	1.4748008395138916	0.0	0.0	5765153.281197977;	%1 onwind
	2	0.0	0.0	3	0.0	1.1963861183763125	0.0	0.0	2136124.274628267;	%1 solar
	2	466996.62931075314	186798.65172430125	3	0.0	2382.6358638303727	0.0	0.0	5941663.7575891465;	%10 CCGT
	2	677412.446313198	270964.9785252792	3	0.0	3456.185950577541	0.0	0.0	2705726.4909073496;	%10 OCGT
	2	4568.621579280366	1827.4486317121462	3	0.0	894.230099682984	0.0	0.0	16684983.921362678;	%10 biomass
	2	0.0	0.0	3	0.0	1.4815853124841887	0.0	0.0	11221801.308137767;	%10 offwind
	2	0.0	0.0	3	0.0	1.4879825561717257	0.0	0.0	5765153.281197978;	%10 onwind
	2	0.0	0.0	3	0.0	1.1893366805473211	0.0	0.0	2136124.274628267;	%10 solar
	2	435754.8075282026	174301.92301128103	3	0.0	2223.238813919401	0.0	0.0	5941663.7575891465;	%11 CCGT
	2	646878.2928675615	258751.31714702462	3	0.0	3300.3994534059257	0.0	0.0	2705726.4909073496;	%11 OCGT
	2	4568.72261860028	1827.4890474401118	3	0.0	894.2498764142258	0.0	0.0	16684983.921362674;	%11 biomass
	2	0.0	0.0	3	0.0	1.518853253905806	0.0	0.0	11157564.772263993;	%11 offwind
	2	0.0	0.0	3	0.0	1.4464086178648843	0.0	0.0	5765153.281197978;	%11 onwind
	2	0.0	0.0	3	0.0	1.20572415356154	0.0	0.0	2136124.274628267;	%11 solar
	2	643392.2178909816	257356.88715639265	3	0.0	3282.6133565866403	0.0	0.0	2705726.4909073496;	%12 OCGT
	2	4568.635787062389	1827.4543148249556	3	0.0	894.2328806150692	0.0	0.0	16684983.921362678;	%12 biomass
	2	0.0	0.0	3	0.0	1.5009571747909594	0.0	0.0	5765153.281197978;	%12 onwind
	2	0.0	0.0	3	0.0	1.2006064345563008	0.0	0.0	2136124.274628267;	%12 solar
	2	483711.052532589	193484.4210130356	3	0.0	2467.9135333295353	0.0	0.0	5941663.7575891465;	%13 CCGT
	2	672863.3364980539	269145.33459922153	3	0.0	3432.976206622724	0.0	0.0	2705726.4909073496;	%13 OCGT
	2	4568.593849770239	1827.437539908096	3	0.0	894.2246721022195	0.0	0.0	16684983.921362678;	%13 biomass
	2	0.0	0.0	3	0.0	1.5257543246987597	0.0	0.0	5765153.281197978;	%13 onwind
	2	0.0	0.0	3	0.0	1.2146185301400887	0.0	0.0	2136124.274628267;	%13 solar
	2	498461.7601163811	199384.70404655242	3	0.0	2543.1722454917403	0.0	0.0	5941663.7575891465;	%14 CCGT
	2	672381.6078735192	268952.6431494077	3	0.0	3430.518407517955	0.0	0.0	2705726.4909073496;	%14 OCGT
	2	4568.406589654658	1827.362635861863	3	0.0	894.188019114241	0.0	0.0	16684983.921362678;	%14 biomass
	2	0.0	0.0	3	0.0	1.4834826213269017	0.0	0.0	11750893.3571214;	%14 offwind
	2	0.0	0.0	3	0.0	1.4688001501300936	0.0	0.0	5765153.281197978;	%14 onwind
	2	0.0	0.0	3	0.0	1.1741918506711848	0.0	0.0	2136124.274628267;	%14 solar
	2	447662.6102274168	179065.04409096672	3	0.0	2283.992909323555	0.0	0.0	5941663.7575891465;	%15 CCGT
	2	627888.3611624561	251155.34446498245	3	0.0	3203.512046747225	0.0	0.0	2705726.4909073496;	%15 OCGT
	2	4568.645931132289	1827.4583724529155	3	0.0	894.2348661445075	0.0	0.0	16684983.921362678;	%15 biomass
	2	0.0	0.0	3	0.0	1.501918165382518	0.0	0.0	12263229.339208078;	%15 offwind
	2	0.0	0.0	3	0.0	1.486723817677117	0.0	0.0	5765153.281197978;	%15 onwind
	2	0.0	0.0	3	0.0	1.1603654162342527	0.0	0.0	2136124.274628267;	%15 solar
	2	468601.95447805006	187440.78179122	3	0.0	2390.8262983573977	0.0	0.0	5941663.7575891465;	%16 CCGT
	2	648368.1799300296	259347.27197201183	3	0.0	3308.000918010355	0.0	0.0	2705726.49090735;	%16 OCGT
	2	4568.535841738451	1827.4143366953808	3	0.0	894.2133180149641	0.0	0.0	16684983.921362681;	%16 biomass
	2	0.0	0.0	3	0.0	1.4842057915620284	0.0	0.0	5765153.281197978;	%16 onwind
	2	0.0	0.0	3	0.0	1.1935390509673542	0.0	0.0	2136124.2746282676;	%16 solar
	2	481783.89360356645	192713.5574414266	3	0.0	2458.0810898141144	0.0	0.0	5941663.7575891465;	%17 CCGT
	2	632054.932081163	252821.9728324652	3	0.0	3224.7700616385864	0.0	0.0	2705726.4909073496;	%17 OCGT
	2	4568.587433066252	1827.4349732265007	3	0.0	894.2234161413685	0.0	0.0	16684983.921362678;	%17 biomass
	2	0.0	0.0	3	0.0	1.5295268384648437	0.0	0.0	5765153.281197978;	%17 onwind
	2	0.0	0.0	3	0.0	1.2230394879356914	0.0	0.0	2136124.2746282676;	%17 solar
	2	505018.1027734795	202007.24110939182	3	0.0	2576.622973334079	0.0	0.0	5941663.7575891465;	%18 CCGT
	2	636973.9672038424	254789.58688153696	3	0.0	3249.8671796114404	0.0	0.0	2705726.4909073496;	%18 OCGT
	2	4568.670124070017	1827.4680496280068	3	0.0	894.2396015012755	0.0	0.0	16684983.921362678;	%18 biomass
	2	0.0	0.0	3	0.0	1.4945757171572027	0.0	0.0	5765153.281197978;	%18 onwind
	2	0.0	0.0	3	0.0	1.1999660538183128	0.0	0.0	2136124.274628267;	%18 solar
	2	473209.0071625708	189283.6028650283	3	0.0	2414.33166919679	0.0	0.0	5941663.7575891465;	%19 CCGT
	2	663008.587829487	265203.43513179483	3	0.0	3382.696876681056	0.0	0.0	2705726.4909073496;	%19 OCGT
	2	4568.564853020283	1827.4259412081135	3	0.0	894.2189964807758	0.0	0.0	16684983.921362678;	%19 biomass
	2	0.0	0.0	3	0.0	1.5248421516882025	0.0	0.0	5765153.281197977;	%19 onwind
	2	0.0	0.0	3	0.0	1.1927498295208692	0.0	0.0	2136124.274628267;	%19 solar
	2	506614.91591100476	202645.96636440195	3	0.0	2584.7699791377795	0.0	0.0	5941663.7575891465;	%2 CCGT
	2	669554.7884747562	267821.9153899025	3	0.0	3416.0958595650823	0.0	0.0	2705726.49090735;	%2 OCGT
	2	4568.542510054601	1827.4170040218403	3	0.0	894.2146232246233	0.0	0.0	16684983.921362678;	%2 biomass
	2	0.0	0.0	3	0.0	1.5027218669013453	0.0	0.0	5765153.281197978;	%2 onwind
	2	0.0	0.0	3	0.0	1.196581099094072	0.0	0.0	2136124.274628267;	%2 solar
	2	552475.2933237642	220990.1173295057	3	0.0	2818.7514965498176	0.0	0.0	5941663.7575891465;	%3 CCGT
	2	687238.0683458183	274895.2273383273	3	0.0	3506.3166752337665	0.0	0.0	2705726.4909073496;	%3 OCGT
	2	4568.409054829089	1827.3636219316356	3	0.0	894.1885016302778	0.0	0.0	16684983.921362678;	%3 biomass
	2	0.0	0.0	3	0.0	1.4864584258945293	0.0	0.0	12208629.695082229;	%3 offwind
	2	0.0	0.0	3	0.0	1.5292150237955453	0.0	0.0	5765153.281197978;	%3 onwind
	2	0.0	0.0	3	0.0	1.198209143669612	0.0	0.0	2136124.274628267;	%3 solar
	2	481599.9181971398	192639.96727885594	3	0.0	2457.1424397813257	0.0	0.0	5941663.7575891465;	%4 CCGT
	2	658628.9911121721	263451.59644486883	3	0.0	3360.351995470266	0.0	0.0	2705726.4909073496;	%4 OCGT
	2	4568.680400770538	1827.4721603082153	3	0.0	894.2416129909059	0.0	0.0	16684983.921362678;	%4 biomass
	2	0.0	0.0	3	0.0	1.4957884576028364	0.0	0.0	5765153.281197978;	%4 onwind
	2	0.0	0.0	3	0.0	1.1726723192802424	0.0	0.0	2136124.274628267;	%4 solar
	2	494665.88257566566	197866.35303026627	3	0.0	2523.805523345233	0.0	0.0	5941663.7575891465;	%5 CCGT
	2	646761.1951767394	258704.4780706957	3	0.0	3299.8020162078533	0.0	0.0	2705726.4909073496;	%5 OCGT
	2	4568.619615015296	1827.4478460061187	3	0.0	894.2297152114497	0.0	0.0	16684983.921362674;	%5 biomass
	2	0.0	0.0	3	0.0	1.4776036426780945	0.0	0.0	5765153.281197978;	%5 onwind
	2	0.0	0.0	3	0.0	1.2010203274751687	0.0	0.0	2136124.274628267;	%5 solar
	2	485015.04247084627	194006.01698833852	3	0.0	2474.5665432186033	0.0	0.0	5941663.7575891465;	%6 CCGT
	2	669861.7091788135	267944.68367152545	3	0.0	3417.661781524559	0.0	0.0	2705726.49090735;	%6 OCGT
	2	4568.695371939177	1827.4781487756707	3	0.0	894.2445433429588	0.0	0.0	16684983.921362678;	%6 biomass
	2	0.0	0.0	3	0.0	1.5005814783479747	0.0	0.0	5765153.281197978;	%6 onwind
	2	0.0	0.0	3	0.0	1.2251314142371388	0.0	0.0	2136124.274628267;	%6 solar
	2	459724.0430925726	183889.61723702904	3	0.0	2345.5308321049624	0.0	0.0	5941663.7575891465;	%7 CCGT
	2	666997.9105959098	266799.1642383639	3	0.0	3403.050564264846	0.0	0.0	2705726.4909073496;	%7 OCGT
	2	4568.6904794847505	1827.4761917939002	3	0.0	894.243585728078	0.0	0.0	16684983.921362674;	%7 biomass
	2	0.0	0.0	3	0.0	1.516994732248627	0.0	0.0	5765153.281197978;	%7 onwind
	2	0.0	0.0	3	0.0	1.197006245263798	0.0	0.0	2136124.2746282676;	%7 solar
	2	519783.0056994758	207913.2022797903	3	0.0	2651.954110711611	0.0	0.0	5941663.7575891465;	%8 CCGT
	2	681633.2563253108	272653.3025301243	3	0.0	3477.7206955373	0.0	0.0	2705726.4909073496;	%8 OCGT
	2	4568.699164274518	1827.479665709807	3	0.0	894.2452856282085	0.0	0.0	16684983.921362678;	%8 biomass
	2	0.0	0.0	3	0.0	1.5334604845819493	0.0	0.0	5765153.281197978;	%8 onwind
	2	0.0	0.0	3	0.0	1.1954783330538625	0.0	0.0	2136124.274628267;	%8 solar
	2	498427.83366042696	199371.1334641708	3	0.0	2542.9991513287086	0.0	0.0	5941663.7575891465;	%9 CCGT
	2	672516.384627307	269006.5538509228	3	0.0	3431.206044016872	0.0	0.0	2705726.4909073496;	%9 OCGT
	2	4568.519863852464	1827.407945540986	3	0.0	894.2101906150841	0.0	0.0	16684983.921362678;	%9 biomass
	2	0.0	0.0	3	0.0	1.4913439307215455	0.0	0.0	5765153.281197977;	%9 onwind
	2	0.0	0.0	3	0.0	1.209411350894525	0.0	0.0	2136124.274628267;	%9 solar
	2	0.0	0.0	3	0.0	0.580085215592837	0.0	0.0	10640712.997169184;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.43506391169462777	0.0	0.0	10640712.997169184;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.631991783365797	0.0	0.0	10640712.997169184;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.47399383752434776	0.0	0.0	10640712.997169184;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.659449536385789	0.0	0.0	10640712.997169184;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.4945871522893417	0.0	0.0	10640712.997169184;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5507066457318421	0.0	0.0	10640712.997169184;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.4130299842988816	0.0	0.0	10640712.997169184;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6551812019467408	0.0	0.0	10640712.997169184;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.49138590146005556	0.0	0.0	10640712.997169184;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6302847597416533	0.0	0.0	10640712.997169184;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.47271356980624	0.0	0.0	10640712.997169184;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5999474269610148	0.0	0.0	10640712.997169184;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.4499605702207611	0.0	0.0	10640712.997169184;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6089915989290302	0.0	0.0	10640712.997169184;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.6021371950905983	0.0	0.0	10640712.997169184;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.4516028963179487	0.0	0.0	10640712.997169184;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5863461194016452	0.0	0.0	10640712.997169184;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.43975958955123384	0.0	0.0	10640712.997169184;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5755688891920677	0.0	0.0	10640712.997169184;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.4316766668940508	0.0	0.0	10640712.997169184;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5511915174265158	0.0	0.0	10640712.997169184;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.4133936380698869	0.0	0.0	10640712.997169184;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6079471250511593	0.0	0.0	10640712.997169184;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.45596034378836947	0.0	0.0	10640712.997169184;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5580657554222411	0.0	0.0	10640712.997169184;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.4185493165666808	0.0	0.0	10640712.997169184;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5614231934911552	0.0	0.0	10640712.997169184;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.42106739511836644	0.0	0.0	10640712.997169184;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.650618105305529	0.0	0.0	10640712.997169184;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.4879635789791468	0.0	0.0	10640712.997169184;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.595590378137482	0.0	0.0	10640712.997169184;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.4466927836031115	0.0	0.0	10640712.997169184;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.619001554918321	0.0	0.0	10640712.997169184;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.46425116618874074	0.0	0.0	10640712.997169184;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6275425364944786	0.0	0.0	10640712.997169184;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.4706569023708589	0.0	0.0	10640712.997169184;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5571835985450139	0.0	0.0	10640712.997169184;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.6485557943224541	0.0	0.0	10640712.997169184;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.48641684574184063	0.0	0.0	10640712.997169184;	%DE1 76 PHS (pump mode)
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

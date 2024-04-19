function mpc = mpc_005_002
%MPC_005_002	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 5 	Weight: 63
%	Time step: 2

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 63;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	-62.541609812363866	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	3721.6923157988067	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	604.1265931787557	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-689.6054933121214	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	832.6238094985399	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	4498.392831935931	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	163.14957176731727	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	1933.1684087182	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	759.385270696113	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	1235.6598560890934	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1143.8555491242214	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	3735.996977172967	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	3912.4580628461777	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2864.5282155093955	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	1334.1909467258242	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	-2505.1898495865576	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	4953.378975232768	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	3265.746167655411	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	2798.9411248941915	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	2846.041230361392	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8.0	8.0	469.0	10	0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.0647722893648326	0.0	0.0	1.0	1.0	276.0	17	1;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.057475217987765	0.0	0.0	43.0	43.0	124.0	38	2;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	23.0	23.0	168.0	28	3;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5.0	5.0	53.0	224	4;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2.0	2.0	53.0	224	5;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.71684375	2.0436	2.8	24.0	24.0	1172.0	10	0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3607083329070386	0.0	0.0	57.0	57.0	309.0	38	2;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	75.0	75.0	419.0	28	3;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2.0	2.0	74.0	224	4;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1.0	1.0	74.0	224	5;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50.0	50.0	1637.0	10	0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.710604488593659	0.0	0.0	20.0	20.0	963.0	17	1;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.853585562070872	0.0	0.0	118.0	118.0	431.0	38	2;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	86.0	86.0	585.0	28	3;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	46.0	224	4;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1.0	1.0	46.0	224	5;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35.0	35.0	1023.0	10	0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.088833939876102	0.0	0.0	7.0	7.0	602.0	17	1;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.767013518770404	0.0	0.0	86.0	86.0	270.0	38	2;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	40.0	40.0	366.0	28	3;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	58.0	224	5;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69.0	69.0	1280.0	10	0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.8471939302646945	0.0	0.0	29.0	29.0	337.0	38	2;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	166.0	166.0	457.0	28	3;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5.0	5.0	72.0	224	4;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3.0	3.0	72.0	224	5;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44.0	44.0	1597.0	10	0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.143788811544169	0.0	0.0	45.0	45.0	421.0	38	2;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	154.0	154.0	571.0	28	3;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	28.0	224	4;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	28.0	224	5;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11.0	11.0	621.0	10	0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.112855845467372	0.0	0.0	34.0	34.0	365.0	17	1;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.35324302848168	0.0	0.0	37.0	37.0	164.0	38	2;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33.0	33.0	222.0	28	3;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2.0	2.0	75.0	224	4;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	75.0	224	5;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54.0	54.0	1672.0	10	0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.312893978598721	0.0	0.0	46.0	46.0	984.0	17	1;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.264427029140186	0.0	0.0	137.0	137.0	440.0	38	2;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56.0	56.0	597.0	28	3;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2.0	2.0	58.0	224	4;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1.0	1.0	58.0	224	5;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52.0	52.0	1279.0	10	0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.6262810248135184	0.0	0.0	18.0	18.0	337.0	38	2;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	173.0	173.0	457.0	28	3;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1.0	1.0	45.0	224	4;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1.0	1.0	45.0	224	5;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12.0	12.0	991.0	10	0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.268137106900414	0.0	0.0	72.0	72.0	261.0	38	2;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	67.0	67.0	354.0	28	3;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2.0	2.0	44.0	224	4;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	44.0	224	5;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18.0	18.0	974.0	10	0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.891844160381921	0.0	0.0	53.0	53.0	257.0	38	2;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	81.0	81.0	348.0	28	3;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17.0	17.0	70.0	224	4;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4.0	4.0	70.0	224	5;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13.0	13.0	1562.0	10	0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7659410193472853	0.0	0.0	49.0	49.0	411.0	38	2;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56.0	56.0	558.0	28	3;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4.0	4.0	81.0	224	4;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1.0	1.0	81.0	224	5;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63.0	63.0	1794.0	10	0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7659792724738406	0.0	0.0	82.0	82.0	472.0	38	2;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	155.0	155.0	641.0	28	3;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5.0	5.0	139.0	224	4;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	139.0	224	5;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65.0	65.0	3109.0	10	0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.82754150730257	0.0	0.0	350.0	350.0	1829.0	17	1;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.008872818827554	0.0	0.0	156.0	156.0	818.0	38	2;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	73.0	73.0	1111.0	28	3;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12.0	12.0	58.0	224	4;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2.0	2.0	58.0	224	5;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38.0	38.0	1279.0	10	0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8097123975580596	0.0	0.0	16.0	16.0	337.0	38	2;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	98.0	98.0	457.0	28	3;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2.0	2.0	35.0	224	4;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1.0	1.0	35.0	224	5;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27.0	27.0	779.0	10	0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.0021012488151015	0.0	0.0	17.0	17.0	205.0	38	2;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	90.0	90.0	278.0	28	3;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11.0	11.0	73.0	224	4;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1.0	1.0	73.0	224	5;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45.0	45.0	1622.0	10	0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7056117115246194	0.0	0.0	70.0	70.0	427.0	38	2;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	91.0	91.0	580.0	28	3;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2.0	2.0	76.0	224	4;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2.0	2.0	76.0	224	5;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54.0	54.0	1682.0	10	0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.1729099679769597	0.0	0.0	109.0	109.0	443.0	38	2;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	111.0	111.0	601.0	28	3;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5.0	5.0	101.0	224	4;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3.0	3.0	101.0	224	5;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60.0	60.0	2241.0	10	0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.420875939036277	0.0	0.0	113.0	113.0	590.0	38	2;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170.0	170.0	801.0	28	3;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14.0	14.0	112.0	224	4;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1.0	1.0	112.0	224	5;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70.0	70.0	2491.0	10	0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.4628304009303057	0.0	0.0	133.0	133.0	656.0	38	2;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	132.0	132.0	890.0	28	3;	%9 solar
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
	2	4797.156050230686	1918.8624200922743	3	0.0	938.9618418928725	0.0	0.0	17519233.117430814;	%0 biomass
	2	0.0	0.0	3	0.0	1.5543928469920079	0.0	0.0	11935315.15781011;	%0 offwind
	2	0.0	0.0	3	0.0	1.5181413803946453	0.0	0.0	6053410.945257876;	%0 onwind
	2	0.0	0.0	3	0.0	1.2007082245169693	0.0	0.0	2242930.4883596804;	%0 solar
	2	559622.8759110206	223849.15036440818	3	0.0	2855.218754648063	0.0	0.0	6238746.945468604;	%1 CCGT
	2	683108.8532865248	273243.5413146099	3	0.0	3485.2492514618607	0.0	0.0	2841012.8154527172;	%1 OCGT
	2	4796.964413419286	1918.7857653677145	3	0.0	938.9243322410033	0.0	0.0	17519233.11743081;	%1 biomass
	2	0.0	0.0	3	0.0	1.5485408814895862	0.0	0.0	6053410.945257876;	%1 onwind
	2	0.0	0.0	3	0.0	1.2562054242951282	0.0	0.0	2242930.4883596804;	%1 solar
	2	490346.4607762908	196138.5843105163	3	0.0	2501.7676570218914	0.0	0.0	6238746.945468604;	%10 CCGT
	2	711283.0686288578	284513.22745154315	3	0.0	3628.995248106418	0.0	0.0	2841012.8154527172;	%10 OCGT
	2	4797.052658244384	1918.8210632977534	3	0.0	938.9416046671332	0.0	0.0	17519233.11743081;	%10 biomass
	2	0.0	0.0	3	0.0	1.555664578108398	0.0	0.0	11782891.373544654;	%10 offwind
	2	0.0	0.0	3	0.0	1.5623816839803122	0.0	0.0	6053410.945257878;	%10 onwind
	2	0.0	0.0	3	0.0	1.2488035145746872	0.0	0.0	2242930.4883596804;	%10 solar
	2	457542.54790461273	183017.0191618451	3	0.0	2334.400754615371	0.0	0.0	6238746.945468604;	%11 CCGT
	2	679222.2075109396	271688.88300437585	3	0.0	3465.419426076222	0.0	0.0	2841012.8154527172;	%11 OCGT
	2	4797.158749530294	1918.8634998121174	3	0.0	938.962370234937	0.0	0.0	17519233.117430806;	%11 biomass
	2	0.0	0.0	3	0.0	1.5947959166010963	0.0	0.0	11715443.010877192;	%11 offwind
	2	0.0	0.0	3	0.0	1.5187290487581286	0.0	0.0	6053410.945257877;	%11 onwind
	2	0.0	0.0	3	0.0	1.266010361239617	0.0	0.0	2242930.4883596804;	%11 solar
	2	675561.8287855308	270224.73151421227	3	0.0	3446.7440244159725	0.0	0.0	2841012.8154527172;	%12 OCGT
	2	4797.067576415508	1918.8270305662033	3	0.0	938.9445246458226	0.0	0.0	17519233.11743081;	%12 biomass
	2	0.0	0.0	3	0.0	1.5760050335305074	0.0	0.0	6053410.945257877;	%12 onwind
	2	0.0	0.0	3	0.0	1.2606367562841156	0.0	0.0	2242930.4883596804;	%12 solar
	2	507896.6051592185	203158.64206368738	3	0.0	2591.309209996012	0.0	0.0	6238746.945468604;	%13 CCGT
	2	706506.5033229566	282602.6013291826	3	0.0	3604.62501695386	0.0	0.0	2841012.8154527172;	%13 OCGT
	2	4797.023542258751	1918.8094169035005	3	0.0	938.9359057073304	0.0	0.0	17519233.11743081;	%13 biomass
	2	0.0	0.0	3	0.0	1.6020420409336977	0.0	0.0	6053410.945257877;	%13 onwind
	2	0.0	0.0	3	0.0	1.275349456647093	0.0	0.0	2242930.4883596804;	%13 solar
	2	523384.8481222001	209353.93924888005	3	0.0	2670.330857766327	0.0	0.0	6238746.945468604;	%14 CCGT
	2	706000.6882671951	282400.27530687803	3	0.0	3602.044327893853	0.0	0.0	2841012.8154527172;	%14 OCGT
	2	4796.826919137391	1918.7307676549563	3	0.0	938.8974200699531	0.0	0.0	17519233.11743081;	%14 biomass
	2	0.0	0.0	3	0.0	1.5576567523932467	0.0	0.0	12338438.02497747;	%14 offwind
	2	0.0	0.0	3	0.0	1.542240157636598	0.0	0.0	6053410.945257878;	%14 onwind
	2	0.0	0.0	3	0.0	1.232901443204744	0.0	0.0	2242930.4883596804;	%14 solar
	2	470045.74073878763	188018.29629551506	3	0.0	2398.1925547897326	0.0	0.0	6238746.945468604;	%15 CCGT
	2	659282.7792205788	263713.1116882316	3	0.0	3363.687649084586	0.0	0.0	2841012.8154527172;	%15 OCGT
	2	4797.078227688903	1918.8312910755612	3	0.0	938.9466094517329	0.0	0.0	17519233.11743081;	%15 biomass
	2	0.0	0.0	3	0.0	1.577014073651644	0.0	0.0	12876390.806168482;	%15 offwind
	2	0.0	0.0	3	0.0	1.561060008560973	0.0	0.0	6053410.945257878;	%15 onwind
	2	0.0	0.0	3	0.0	1.2183836870459654	0.0	0.0	2242930.4883596804;	%15 solar
	2	492032.05220195255	196812.82088078102	3	0.0	2510.367613275268	0.0	0.0	6238746.945468604;	%16 CCGT
	2	680786.588926531	272314.6355706124	3	0.0	3473.400963910872	0.0	0.0	2841012.8154527177;	%16 OCGT
	2	4796.962633825375	1918.7850535301498	3	0.0	938.9239839157123	0.0	0.0	17519233.117430814;	%16 biomass
	2	0.0	0.0	3	0.0	1.5584160811401297	0.0	0.0	6053410.945257877;	%16 onwind
	2	0.0	0.0	3	0.0	1.253216003515722	0.0	0.0	2242930.488359681;	%16 solar
	2	505873.0882837448	202349.23531349792	3	0.0	2580.98514430482	0.0	0.0	6238746.945468604;	%17 CCGT
	2	663657.6786852211	265463.0714740885	3	0.0	3386.0085647205156	0.0	0.0	2841012.8154527172;	%17 OCGT
	2	4797.016804719564	1918.8067218878257	3	0.0	938.9345869484368	0.0	0.0	17519233.11743081;	%17 biomass
	2	0.0	0.0	3	0.0	1.606003180388086	0.0	0.0	6053410.945257878;	%17 onwind
	2	0.0	0.0	3	0.0	1.284191462332476	0.0	0.0	2242930.488359681;	%17 solar
	2	530269.0079121535	212107.6031648614	3	0.0	2705.4541220007827	0.0	0.0	6238746.945468604;	%18 CCGT
	2	668822.6655640345	267529.0662256138	3	0.0	3412.3605385920127	0.0	0.0	2841012.8154527172;	%18 OCGT
	2	4797.103630273517	1918.841452109407	3	0.0	938.9515815763393	0.0	0.0	17519233.11743081;	%18 biomass
	2	0.0	0.0	3	0.0	1.5693045030150627	0.0	0.0	6053410.945257877;	%18 onwind
	2	0.0	0.0	3	0.0	1.2599643565092284	0.0	0.0	2242930.4883596804;	%18 solar
	2	496869.45752069936	198747.78300827974	3	0.0	2535.0482526566293	0.0	0.0	6238746.945468604;	%19 CCGT
	2	696159.0172209614	278463.60688838456	3	0.0	3551.831720515109	0.0	0.0	2841012.8154527172;	%19 OCGT
	2	4796.993095671298	1918.797238268519	3	0.0	938.9299463048145	0.0	0.0	17519233.11743081;	%19 biomass
	2	0.0	0.0	3	0.0	1.6010842592726127	0.0	0.0	6053410.945257876;	%19 onwind
	2	0.0	0.0	3	0.0	1.2523873209969127	0.0	0.0	2242930.4883596804;	%19 solar
	2	531945.661706555	212778.26468262204	3	0.0	2714.0084780946686	0.0	0.0	6238746.945468604;	%2 CCGT
	2	703032.5278984939	281213.0111593976	3	0.0	3586.900652543336	0.0	0.0	2841012.8154527177;	%2 OCGT
	2	4796.969635557331	1918.7878542229323	3	0.0	938.9253543858546	0.0	0.0	17519233.11743081;	%2 biomass
	2	0.0	0.0	3	0.0	1.5778579602464127	0.0	0.0	6053410.945257877;	%2 onwind
	2	0.0	0.0	3	0.0	1.2564101540487755	0.0	0.0	2242930.4883596804;	%2 solar
	2	580099.0579899525	232039.623195981	3	0.0	2959.6890713773087	0.0	0.0	6238746.945468604;	%3 CCGT
	2	721599.9717631093	288639.9887052437	3	0.0	3681.6325089954553	0.0	0.0	2841012.8154527172;	%3 OCGT
	2	4796.829507570544	1918.7318030282174	3	0.0	938.8979267117917	0.0	0.0	17519233.11743081;	%3 biomass
	2	0.0	0.0	3	0.0	1.5607813471892555	0.0	0.0	12819061.17983634;	%3 offwind
	2	0.0	0.0	3	0.0	1.6056757749853225	0.0	0.0	6053410.945257878;	%3 onwind
	2	0.0	0.0	3	0.0	1.2581196008530926	0.0	0.0	2242930.4883596804;	%3 solar
	2	505679.91410699685	202271.96564279872	3	0.0	2579.9995617703917	0.0	0.0	6238746.945468604;	%4 CCGT
	2	691560.4406677807	276624.1762671123	3	0.0	3528.3695952437793	0.0	0.0	2841012.8154527172;	%4 OCGT
	2	4797.114420809065	1918.8457683236259	3	0.0	938.9536936404511	0.0	0.0	17519233.11743081;	%4 biomass
	2	0.0	0.0	3	0.0	1.5705778804829782	0.0	0.0	6053410.945257877;	%4 onwind
	2	0.0	0.0	3	0.0	1.2313059352442546	0.0	0.0	2242930.4883596804;	%4 solar
	2	519399.17670444894	207759.67068177956	3	0.0	2649.9957995124946	0.0	0.0	6238746.945468604;	%5 CCGT
	2	679099.2549355762	271639.7019742305	3	0.0	3464.792117018246	0.0	0.0	2841012.8154527172;	%5 OCGT
	2	4797.050595766061	1918.8202383064245	3	0.0	938.9412009720222	0.0	0.0	17519233.117430806;	%5 biomass
	2	0.0	0.0	3	0.0	1.551483824811999	0.0	0.0	6053410.945257877;	%5 onwind
	2	0.0	0.0	3	0.0	1.2610713438489272	0.0	0.0	2242930.4883596804;	%5 solar
	2	509265.7945943886	203706.31783775543	3	0.0	2598.2948703795337	0.0	0.0	6238746.945468604;	%6 CCGT
	2	703354.7946377542	281341.9178551017	3	0.0	3588.544870600787	0.0	0.0	2841012.8154527177;	%6 OCGT
	2	4797.1301405361355	1918.8520562144543	3	0.0	938.9567705101067	0.0	0.0	17519233.11743081;	%6 biomass
	2	0.0	0.0	3	0.0	1.5756105522653734	0.0	0.0	6053410.945257877;	%6 onwind
	2	0.0	0.0	3	0.0	1.2863879849489956	0.0	0.0	2242930.4883596804;	%6 solar
	2	482710.24524720124	193084.0980988805	3	0.0	2462.8073737102104	0.0	0.0	6238746.945468604;	%7 CCGT
	2	700347.8061257054	280139.1224502821	3	0.0	3573.203092478088	0.0	0.0	2841012.8154527172;	%7 OCGT
	2	4797.1250034589875	1918.8500013835953	3	0.0	938.9557650144819	0.0	0.0	17519233.117430806;	%7 biomass
	2	0.0	0.0	3	0.0	1.5928444688610583	0.0	0.0	6053410.945257877;	%7 onwind
	2	0.0	0.0	3	0.0	1.2568565575269879	0.0	0.0	2242930.488359681;	%7 solar
	2	545772.1559844497	218308.86239377983	3	0.0	2784.5518162471913	0.0	0.0	6238746.945468604;	%8 CCGT
	2	715714.9191415764	286285.96765663056	3	0.0	3651.606730314165	0.0	0.0	2841012.8154527172;	%8 OCGT
	2	4797.1341224882435	1918.8536489952974	3	0.0	938.9575499096189	0.0	0.0	17519233.11743081;	%8 biomass
	2	0.0	0.0	3	0.0	1.6101335088110469	0.0	0.0	6053410.945257877;	%8 onwind
	2	0.0	0.0	3	0.0	1.2552522497065555	0.0	0.0	2242930.4883596804;	%8 solar
	2	523349.2253434483	209339.6901373793	3	0.0	2670.1491088951443	0.0	0.0	6238746.945468604;	%9 CCGT
	2	706142.2038586723	282456.8815434689	3	0.0	3602.7663462177156	0.0	0.0	2841012.8154527172;	%9 OCGT
	2	4796.945857045088	1918.7783428180353	3	0.0	938.9207001458383	0.0	0.0	17519233.11743081;	%9 biomass
	2	0.0	0.0	3	0.0	1.5659111272576227	0.0	0.0	6053410.945257876;	%9 onwind
	2	0.0	0.0	3	0.0	1.2698819184392511	0.0	0.0	2242930.4883596804;	%9 solar
	2	0.0	0.0	3	0.0	0.6090894763724788	0.0	0.0	11172748.647027643;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.4568171072793591	0.0	0.0	11172748.647027643;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6635913725340868	0.0	0.0	11172748.647027643;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.4976935294005651	0.0	0.0	11172748.647027643;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6924220132050783	0.0	0.0	11172748.647027643;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.5193165099038087	0.0	0.0	11172748.647027643;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5782419780184342	0.0	0.0	11172748.647027643;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.4336814835138257	0.0	0.0	11172748.647027643;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6879402620440778	0.0	0.0	11172748.647027643;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.5159551965330583	0.0	0.0	11172748.647027643;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6617989977287361	0.0	0.0	11172748.647027643;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.496349248296552	0.0	0.0	11172748.647027643;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6299447983090656	0.0	0.0	11172748.647027643;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.4724585987317992	0.0	0.0	11172748.647027643;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6394411788754818	0.0	0.0	11172748.647027643;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.6322440548451282	0.0	0.0	11172748.647027643;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.4741830411338462	0.0	0.0	11172748.647027643;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6156634253717274	0.0	0.0	11172748.647027643;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.46174756902879555	0.0	0.0	11172748.647027643;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.604347333651671	0.0	0.0	11172748.647027643;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.4532605002387533	0.0	0.0	11172748.647027643;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5787510932978416	0.0	0.0	11172748.647027643;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.4340633199733812	0.0	0.0	11172748.647027643;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6383444813037172	0.0	0.0	11172748.647027643;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.4787583609777879	0.0	0.0	11172748.647027643;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5859690431933532	0.0	0.0	11172748.647027643;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.43947678239501486	0.0	0.0	11172748.647027643;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.589494353165713	0.0	0.0	11172748.647027643;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.44212076487428476	0.0	0.0	11172748.647027643;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6831490105708055	0.0	0.0	11172748.647027643;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.5123617579281041	0.0	0.0	11172748.647027643;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.625369897044356	0.0	0.0	11172748.647027643;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.469027422783267	0.0	0.0	11172748.647027643;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6499516326642371	0.0	0.0	11172748.647027643;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.48746372449817776	0.0	0.0	11172748.647027643;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.6589196633192025	0.0	0.0	11172748.647027643;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.49418974748940186	0.0	0.0	11172748.647027643;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5850427784722646	0.0	0.0	11172748.647027643;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.6809835840385768	0.0	0.0	11172748.647027643;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.5107376880289326	0.0	0.0	11172748.647027643;	%DE1 76 PHS (pump mode)
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

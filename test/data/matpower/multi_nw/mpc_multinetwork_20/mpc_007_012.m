function mpc = mpc_007_012
%MPC_007_012	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 7 	Weight: 39
%	Time step: 12

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1672.35990413242	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	5578.85855308611	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1888.5728758864086	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	2478.4078080536206	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	2493.0020496671127	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	6755.507649707436	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	337.4094199653608	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2920.940511730517	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	3821.650826997855	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	3673.292771330037	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	2648.4043647326657	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	3060.5313342133004	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5761.580574791317	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2972.3766518345956	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3746.663173964089	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	3301.3705085947454	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	3077.547023891747	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5193.751132747236	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5488.1055792322095	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4244.751863614404	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8.0	8.0	469.0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.11583222921456	0.0	0.0	1.0	1.0	276.0;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.19874311950359	0.0	0.0	43.0	43.0	124.0;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.018682889800027	0.0	0.0	23.0	23.0	168.0;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5.0	5.0	53.0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2.0	2.0	53.0;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.716843749999999	2.0436	2.8	24.0	24.0	1172.0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.725869379940092	0.0	0.0	57.0	57.0	309.0;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.699709726720505	0.0	0.0	75.0	75.0	419.0;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2.0	2.0	74.0;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1.0	1.0	74.0;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999999	2.0436	2.8	50.0	50.0	1637.0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.743195391260197	0.0	0.0	20.0	20.0	963.0;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.899485058892019	0.0	0.0	118.0	118.0	431.0;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.080350888210747	0.0	0.0	86.0	86.0	585.0;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	46.0;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1.0	1.0	46.0;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35.0	35.0	1023.0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.6992816950928775	0.0	0.0	7.0	7.0	602.0;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.286032183864318	0.0	0.0	86.0	86.0	270.0;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.475094212898853	0.0	0.0	40.0	40.0	366.0;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	58.0;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69.0	69.0	1280.0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8560623853500497	0.0	0.0	29.0	29.0	337.0;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.741669863549268	0.0	0.0	166.0	166.0	457.0;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.292	78.4	5.0	5.0	5.0	72.0;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3.0	3.0	72.0;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44.0	44.0	1597.0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.403043535055404	0.0	0.0	45.0	45.0	421.0;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.1875351740477	0.0	0.0	154.0	154.0	571.0;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	28.0;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	28.0;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11.0	11.0	621.0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.045578534058372	0.0	0.0	34.0	34.0	365.0;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.220232603286098	0.0	0.0	37.0	37.0	164.0;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.086331087339863	0.0	0.0	33.0	33.0	222.0;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2.0	2.0	75.0;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	75.0;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54.0	54.0	1672.0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.53335839827153	0.0	0.0	46.0	46.0	984.0;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.812487301209975	0.0	0.0	137.0	137.0	440.0;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.389298195742072	0.0	0.0	56.0	56.0	597.0;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2.0	2.0	58.0;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1.0	1.0	58.0;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52.0	52.0	1279.0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.066905340154448	0.0	0.0	18.0	18.0	337.0;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.486189741737551	0.0	0.0	173.0	173.0	457.0;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1.0	1.0	45.0;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1.0	1.0	45.0;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333334	2.0436	2.8	12.0	12.0	991.0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.18383993027638	0.0	0.0	72.0	72.0	261.0;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.966134093311785	0.0	0.0	67.0	67.0	354.0;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2.0	2.0	44.0;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	44.0;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18.0	18.0	974.0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.348256726759682	0.0	0.0	53.0	53.0	257.0;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.525052962906769	0.0	0.0	81.0	81.0	348.0;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17.0	17.0	70.0;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4.0	4.0	70.0;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13.0	13.0	1562.0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.086446197633803	0.0	0.0	49.0	49.0	411.0;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.92910685004568	0.0	0.0	56.0	56.0	558.0;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4.0	4.0	81.0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1.0	1.0	81.0;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63.0	63.0	1794.0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.625640278910398	0.0	0.0	82.0	82.0	472.0;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.466269106679217	0.0	0.0	155.0	155.0	641.0;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999999	78.4	5.0	5.0	5.0	139.0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	139.0;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65.0	65.0	3109.0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.30036023626575	0.0	0.0	350.0	350.0	1829.0;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.239682378772095	0.0	0.0	156.0	156.0	818.0;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.329997949424419	0.0	0.0	73.0	73.0	1111.0;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12.0	12.0	58.0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2.0	2.0	58.0;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38.0	38.0	1279.0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.1580835350116567	0.0	0.0	16.0	16.0	337.0;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.9230774586666	0.0	0.0	98.0	98.0	457.0;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2.0	2.0	35.0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1.0	1.0	35.0;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27.0	27.0	779.0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.9347619267394807	0.0	0.0	17.0	17.0	205.0;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.844362694631712	0.0	0.0	90.0	90.0	278.0;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11.0	11.0	73.0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1.0	1.0	73.0;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45.0	45.0	1622.0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.968319985087527	0.0	0.0	70.0	70.0	427.0;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.65772733482649	0.0	0.0	91.0	91.0	580.0;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2.0	2.0	76.0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26	78.4	2.8	2.0	2.0	76.0;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54.0	54.0	1682.0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.202220676723227	0.0	0.0	109.0	109.0	443.0;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.330942183444586	0.0	0.0	111.0	111.0	601.0;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5.0	5.0	101.0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3.0	3.0	101.0;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60.0	60.0	2241.0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.219070325282285	0.0	0.0	113.0	113.0	590.0;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.7669012128533	0.0	0.0	170.0	170.0	801.0;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14.0	14.0	112.0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1.0	1.0	112.0;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70.0	70.0	2491.0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.625192278909553	0.0	0.0	133.0	133.0	656.0;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.531490039186341	0.0	0.0	132.0	132.0	890.0;	%9 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1776.8	266.52	3.5	1	1	1;	%DE1 0 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-266.52	-1776.8	3.5	1	1	1;	%DE1 0 PHS (pump mode)
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	1294.0	194.1	3.5	1	1	1;	%DE1 12 PHS
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	-194.1	-1294.0	3.5	1	1	1;	%DE1 12 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.42	1.113	3.5	1	1	1;	%DE1 14 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-1.113	-7.42	3.5	1	1	1;	%DE1 14 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.0	33.0	3.5	1	1	1;	%DE1 17 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-33.0	-220.0	3.5	1	1	1;	%DE1 17 PHS (pump mode)
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	960.0	144.0	3.5	1	1	1;	%DE1 2 PHS
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	-144.0	-960.0	3.5	1	1	1;	%DE1 2 PHS (pump mode)
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	39.8	5.97	3.5	1	1	1;	%DE1 20 PHS
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	-5.97	-39.8	3.5	1	1	1;	%DE1 20 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	160.0	24.0	3.5	1	1	1;	%DE1 21 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.0	-160.0	3.5	1	1	1;	%DE1 21 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	124.0	18.599999999999998	3.5	1	1	1;	%DE1 25 hydro
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	167.1	25.064999999999998	3.5	1	1	1;	%DE1 27 PHS
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	-25.064999999999998	-167.1	3.5	1	1	1;	%DE1 27 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	109.0	16.349999999999998	3.5	1	1	1;	%DE1 31 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-16.349999999999998	-109.0	3.5	1	1	1;	%DE1 31 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	480.0	72.0	3.5	1	1	1;	%DE1 33 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-72.0	-480.0	3.5	1	1	1;	%DE1 33 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	92.0	13.799999999999999	3.5	1	1	1;	%DE1 36 PHS
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.799999999999999	-92.0	3.5	1	1	1;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	278.0	41.699999999999996	3.5	1	1	1;	%DE1 43 PHS
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	-41.699999999999996	-278.0	3.5	1	1	1;	%DE1 43 PHS (pump mode)
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	164.0	24.599999999999998	3.5	1	1	1;	%DE1 49 PHS
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.599999999999998	-164.0	3.5	1	1	1;	%DE1 49 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.7	11.955	3.5	1	1	1;	%DE1 50 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-11.955	-79.7	3.5	1	1	1;	%DE1 50 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	399.8	59.97	3.5	1	1	1;	%DE1 52 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-59.97	-399.8	3.5	1	1	1;	%DE1 52 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	162.0	24.3	3.5	1	1	1;	%DE1 54 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.3	-162.0	3.5	1	1	1;	%DE1 54 PHS (pump mode)
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	90.0	13.5	3.5	1	1	1;	%DE1 61 PHS
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.5	-90.0	3.5	1	1	1;	%DE1 61 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	1060.0	159.0	3.5	1	1	1;	%DE1 7 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-159.0	-1060.0	3.5	1	1	1;	%DE1 7 PHS (pump mode)
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	45.5	6.825	3.5	1	1	1;	%DE1 71 hydro
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	18.0	3.5	1	1	1;	%DE1 76 PHS
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	-18.0	-120.0	3.5	1	1	1;	%DE1 76 PHS (pump mode)
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0	investment
mpc.gencost = [
	2	2969.6680310951865	1187.8672124380746	3	0.0	581.2620926003497	0.0	0.0	10845239.548885742;	%0 biomass
	2	0.0	0.0	3	0.0	0.9622431909950524	0.0	0.0	7388528.431025307;	%0 offwind
	2	0.0	0.0	3	0.0	0.9398018069109709	0.0	0.0	3747349.632778685;	%0 onwind
	2	0.0	0.0	3	0.0	0.743295567558124	0.0	0.0	1388480.7785083735;	%0 solar
	2	346433.20889729843	138573.28355891936	3	0.0	1767.5163719249915	0.0	0.0	3862081.442432945;	%1 CCGT
	2	422876.90917737246	169150.76367094897	3	0.0	2157.5352509049612	0.0	0.0	1758722.2190897774;	%1 OCGT
	2	2969.5493987833675	1187.819759513347	3	0.0	581.2388723396688	0.0	0.0	10845239.54888574;	%1 biomass
	2	0.0	0.0	3	0.0	0.9586205456840295	0.0	0.0	3747349.632778685;	%1 onwind
	2	0.0	0.0	3	0.0	0.7776509769446032	0.0	0.0	1388480.7785083735;	%1 solar
	2	303547.8090519896	121419.12362079581	3	0.0	1548.7133114897424	0.0	0.0	3862081.442432945;	%10 CCGT
	2	440318.0901035787	176127.2360414315	3	0.0	2246.5208678754016	0.0	0.0	1758722.2190897774;	%10 OCGT
	2	2969.6040265322376	1187.841610612895	3	0.0	581.2495647939396	0.0	0.0	10845239.54888574;	%10 biomass
	2	0.0	0.0	3	0.0	0.9630304531147226	0.0	0.0	7294170.850289548;	%10 offwind
	2	0.0	0.0	3	0.0	0.9671886615116218	0.0	0.0	3747349.632778686;	%10 onwind
	2	0.0	0.0	3	0.0	0.7730688423557587	0.0	0.0	1388480.7785083735;	%10 solar
	2	283240.6248933317	113296.24995733268	3	0.0	1445.1052290476107	0.0	0.0	3862081.442432945;	%11 CCGT
	2	420470.89036391495	168188.35614556598	3	0.0	2145.2596447138517	0.0	0.0	1758722.2190897774;	%11 OCGT
	2	2969.669702090182	1187.8678808360726	3	0.0	581.2624196692467	0.0	0.0	10845239.548885738;	%11 biomass
	2	0.0	0.0	3	0.0	0.9872546150387739	0.0	0.0	7252417.101971595;	%11 offwind
	2	0.0	0.0	3	0.0	0.9401656016121748	0.0	0.0	3747349.6327786855;	%11 onwind
	2	0.0	0.0	3	0.0	0.7837206998150009	0.0	0.0	1388480.7785083735;	%11 solar
	2	418204.9416291381	167281.97665165522	3	0.0	2133.6986817813163	0.0	0.0	1758722.2190897774;	%12 OCGT
	2	2969.6132615905526	1187.845304636221	3	0.0	581.251372399795	0.0	0.0	10845239.54888574;	%12 biomass
	2	0.0	0.0	3	0.0	0.9756221636141236	0.0	0.0	3747349.6327786855;	%12 onwind
	2	0.0	0.0	3	0.0	0.7803941824615954	0.0	0.0	1388480.7785083735;	%12 solar
	2	314412.18414618284	125764.87365847314	3	0.0	1604.1437966641981	0.0	0.0	3862081.442432945;	%13 CCGT
	2	437361.168723735	174944.467489494	3	0.0	2231.4345343047703	0.0	0.0	1758722.2190897774;	%13 OCGT
	2	2969.5860023506557	1187.8344009402622	3	0.0	581.2460368664426	0.0	0.0	10845239.54888574;	%13 biomass
	2	0.0	0.0	3	0.0	0.9917403110541938	0.0	0.0	3747349.6327786855;	%13 onwind
	2	0.0	0.0	3	0.0	0.7895020445910577	0.0	0.0	1388480.7785083735;	%13 solar
	2	324000.1440756477	129600.05763025908	3	0.0	1653.0619595696312	0.0	0.0	3862081.442432945;	%14 CCGT
	2	437048.04511778743	174819.218047115	3	0.0	2229.8369648866706	0.0	0.0	1758722.2190897774;	%14 OCGT
	2	2969.4642832755276	1187.785713310211	3	0.0	581.2222124242567	0.0	0.0	10845239.54888574;	%14 biomass
	2	0.0	0.0	3	0.0	0.9642637038624862	0.0	0.0	7638080.68212891;	%14 offwind
	2	0.0	0.0	3	0.0	0.9547200975845608	0.0	0.0	3747349.632778686;	%14 onwind
	2	0.0	0.0	3	0.0	0.7632247029362701	0.0	0.0	1388480.7785083735;	%14 solar
	2	290980.69664782094	116392.27865912837	3	0.0	1484.5953910603107	0.0	0.0	3862081.442432945;	%15 CCGT
	2	408127.43475559645	163250.9739022386	3	0.0	2082.2828303856963	0.0	0.0	1758722.2190897774;	%15 OCGT
	2	2969.6198552359874	1187.847942094395	3	0.0	581.2526629939299	0.0	0.0	10845239.54888574;	%15 biomass
	2	0.0	0.0	3	0.0	0.9762468074986367	0.0	0.0	7971099.07048525;	%15 offwind
	2	0.0	0.0	3	0.0	0.966370481490126	0.0	0.0	3747349.632778686;	%15 onwind
	2	0.0	0.0	3	0.0	0.7542375205522643	0.0	0.0	1388480.7785083735;	%15 solar
	2	304591.27041073254	121836.50816429302	3	0.0	1554.0370939323086	0.0	0.0	3862081.442432945;	%16 CCGT
	2	421439.31695451925	168575.7267818077	3	0.0	2150.2005967067307	0.0	0.0	1758722.2190897777;	%16 OCGT
	2	2969.5482971299934	1187.8193188519974	3	0.0	581.2386567097267	0.0	0.0	10845239.548885742;	%16 biomass
	2	0.0	0.0	3	0.0	0.9647337645153184	0.0	0.0	3747349.6327786855;	%16 onwind
	2	0.0	0.0	3	0.0	0.7758003831287803	0.0	0.0	1388480.7785083738;	%16 solar
	2	313159.5308423182	125263.81233692728	3	0.0	1597.7527083791742	0.0	0.0	3862081.442432945;	%17 CCGT
	2	410835.705852756	164334.28234110237	3	0.0	2096.100540065081	0.0	0.0	1758722.2190897774;	%17 OCGT
	2	2969.5818314930634	1187.8327325972255	3	0.0	581.2452204918895	0.0	0.0	10845239.54888574;	%17 biomass
	2	0.0	0.0	3	0.0	0.9941924450021484	0.0	0.0	3747349.632778686;	%17 onwind
	2	0.0	0.0	3	0.0	0.7949756671581993	0.0	0.0	1388480.7785083738;	%17 solar
	2	328261.76680276165	131304.70672110468	3	0.0	1674.8049326671512	0.0	0.0	3862081.442432945;	%18 CCGT
	2	414033.07868249755	165613.23147299903	3	0.0	2112.4136667474363	0.0	0.0	1758722.2190897774;	%18 OCGT
	2	2969.635580645511	1187.8542322582043	3	0.0	581.2557409758291	0.0	0.0	10845239.54888574;	%18 biomass
	2	0.0	0.0	3	0.0	0.9714742161521818	0.0	0.0	3747349.6327786855;	%18 onwind
	2	0.0	0.0	3	0.0	0.7799779349819033	0.0	0.0	1388480.7785083735;	%18 solar
	2	307585.854655671	123034.34186226841	3	0.0	1569.3155849779132	0.0	0.0	3862081.442432945;	%19 CCGT
	2	430955.5820891666	172382.23283566663	3	0.0	2198.7529698426865	0.0	0.0	1758722.2190897774;	%19 OCGT
	2	2969.567154463184	1187.8268617852736	3	0.0	581.2423477125043	0.0	0.0	10845239.54888574;	%19 biomass
	2	0.0	0.0	3	0.0	0.9911473985973317	0.0	0.0	3747349.632778685;	%19 onwind
	2	0.0	0.0	3	0.0	0.775287389188565	0.0	0.0	1388480.7785083735;	%19 solar
	2	329299.69534215314	131719.87813686125	3	0.0	1680.1004864395568	0.0	0.0	3862081.442432945;	%2 CCGT
	2	435210.6125085915	174084.2450034366	3	0.0	2220.4623087173036	0.0	0.0	1758722.2190897777;	%2 OCGT
	2	2969.5526315354905	1187.821052614196	3	0.0	581.2395050960051	0.0	0.0	10845239.54888574;	%2 biomass
	2	0.0	0.0	3	0.0	0.9767692134858745	0.0	0.0	3747349.6327786855;	%2 onwind
	2	0.0	0.0	3	0.0	0.7777777144111467	0.0	0.0	1388480.7785083735;	%2 solar
	2	359108.94066044677	143643.57626417873	3	0.0	1832.1884727573815	0.0	0.0	3862081.442432945;	%3 CCGT
	2	446704.74442478194	178681.89776991276	3	0.0	2279.1058389019486	0.0	0.0	1758722.2190897774;	%3 OCGT
	2	2969.465885638908	1187.786354255563	3	0.0	581.2225260596806	0.0	0.0	10845239.54888574;	%3 biomass
	2	0.0	0.0	3	0.0	0.9661979768314439	0.0	0.0	7935609.301803448;	%3 offwind
	2	0.0	0.0	3	0.0	0.9939897654671044	0.0	0.0	3747349.632778686;	%3 onwind
	2	0.0	0.0	3	0.0	0.7788359433852478	0.0	0.0	1388480.7785083735;	%3 solar
	2	313039.9468281409	125215.97873125636	3	0.0	1597.1425858578616	0.0	0.0	3862081.442432945;	%4 CCGT
	2	428108.84422291187	171243.53768916475	3	0.0	2184.2287970556727	0.0	0.0	1758722.2190897774;	%4 OCGT
	2	2969.6422605008497	1187.8569042003398	3	0.0	581.2570484440888	0.0	0.0	10845239.54888574;	%4 biomass
	2	0.0	0.0	3	0.0	0.9722624974418437	0.0	0.0	3747349.6327786855;	%4 onwind
	2	0.0	0.0	3	0.0	0.7622370075321576	0.0	0.0	1388480.7785083735;	%4 solar
	2	321532.8236741827	128613.12946967308	3	0.0	1640.4735901744014	0.0	0.0	3862081.442432945;	%5 CCGT
	2	420394.7768648806	168157.91074595222	3	0.0	2144.8713105351044	0.0	0.0	1758722.2190897774;	%5 OCGT
	2	2969.6027497599425	1187.841099903977	3	0.0	581.2493148874423	0.0	0.0	10845239.548885738;	%5 biomass
	2	0.0	0.0	3	0.0	0.9604423677407614	0.0	0.0	3747349.6327786855;	%5 onwind
	2	0.0	0.0	3	0.0	0.7806632128588596	0.0	0.0	1388480.7785083735;	%5 solar
	2	315259.7776060501	126103.91104242003	3	0.0	1608.4682530920923	0.0	0.0	3862081.442432945;	%6 CCGT
	2	435410.1109662288	174164.04438649153	3	0.0	2221.480157990963	0.0	0.0	1758722.2190897777;	%6 OCGT
	2	2969.6519917604646	1187.860796704186	3	0.0	581.2589531729233	0.0	0.0	10845239.54888574;	%6 biomass
	2	0.0	0.0	3	0.0	0.9753779609261836	0.0	0.0	3747349.6327786855;	%6 onwind
	2	0.0	0.0	3	0.0	0.7963354192541402	0.0	0.0	1388480.7785083735;	%6 solar
	2	298820.62801017216	119528.25120406887	3	0.0	1524.5950408682254	0.0	0.0	3862081.442432945;	%7 CCGT
	2	433548.6418873414	173419.45675493655	3	0.0	2211.98286677215	0.0	0.0	1758722.2190897774;	%7 OCGT
	2	2969.6488116650876	1187.8595246660352	3	0.0	581.2583307232507	0.0	0.0	10845239.548885738;	%7 biomass
	2	0.0	0.0	3	0.0	0.9860465759616076	0.0	0.0	3747349.6327786855;	%7 onwind
	2	0.0	0.0	3	0.0	0.7780540594214687	0.0	0.0	1388480.7785083738;	%7 solar
	2	337858.9537046593	135143.58148186372	3	0.0	1723.770171962547	0.0	0.0	3862081.442432945;	%8 CCGT
	2	443061.616611452	177224.64664458082	3	0.0	2260.518452099245	0.0	0.0	1758722.2190897774;	%8 OCGT
	2	2969.6544567784363	1187.8617827113744	3	0.0	581.2594356583355	0.0	0.0	10845239.54888574;	%8 biomass
	2	0.0	0.0	3	0.0	0.996749314978267	0.0	0.0	3747349.6327786855;	%8 onwind
	2	0.0	0.0	3	0.0	0.7770609164850106	0.0	0.0	1388480.7785083735;	%8 solar
	2	323978.0918792775	129591.23675171101	3	0.0	1652.9494483636606	0.0	0.0	3862081.442432945;	%9 CCGT
	2	437135.6500077495	174854.2600030998	3	0.0	2230.2839286109665	0.0	0.0	1758722.2190897774;	%9 OCGT
	2	2969.537911504102	1187.8151646016408	3	0.0	581.2366238998046	0.0	0.0	10845239.54888574;	%9 biomass
	2	0.0	0.0	3	0.0	0.9693735549690046	0.0	0.0	3747349.632778685;	%9 onwind
	2	0.0	0.0	3	0.0	0.7861173780814412	0.0	0.0	1388480.7785083735;	%9 solar
	2	0.0	0.0	3	0.0	0.377055390135344	0.0	0.0	6916463.448159969;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.377055390135344	0.0	0.0	-6916463.448159969;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.41079465918776803	0.0	0.0	6916463.448159969;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.41079465918776803	0.0	0.0	-6916463.448159969;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4286421986507628	0.0	0.0	6916463.448159969;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.4286421986507628	0.0	0.0	-6916463.448159969;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3579593197256974	0.0	0.0	6916463.448159969;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.3579593197256974	0.0	0.0	-6916463.448159969;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4258677812653815	0.0	0.0	6916463.448159969;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.4258677812653815	0.0	0.0	-6916463.448159969;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4096850938320747	0.0	0.0	6916463.448159969;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.4096850938320747	0.0	0.0	-6916463.448159969;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38996582752465964	0.0	0.0	6916463.448159969;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.38996582752465964	0.0	0.0	-6916463.448159969;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39584453930386965	0.0	0.0	6916463.448159969;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.3913891768088889	0.0	0.0	6916463.448159969;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.3913891768088889	0.0	0.0	-6916463.448159969;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38112497761106934	0.0	0.0	6916463.448159969;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.38112497761106934	0.0	0.0	-6916463.448159969;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.374119777974844	0.0	0.0	6916463.448159969;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.374119777974844	0.0	0.0	-6916463.448159969;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3582744863272353	0.0	0.0	6916463.448159969;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.3582744863272353	0.0	0.0	-6916463.448159969;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39516563128325355	0.0	0.0	6916463.448159969;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.39516563128325355	0.0	0.0	-6916463.448159969;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3627427410244567	0.0	0.0	6916463.448159969;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.3627427410244567	0.0	0.0	-6916463.448159969;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3649250757692509	0.0	0.0	6916463.448159969;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.3649250757692509	0.0	0.0	-6916463.448159969;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.42290176844859384	0.0	0.0	6916463.448159969;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.42290176844859384	0.0	0.0	-6916463.448159969;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38713374578936327	0.0	0.0	6916463.448159969;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.38713374578936327	0.0	0.0	-6916463.448159969;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40235101069690865	0.0	0.0	6916463.448159969;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.40235101069690865	0.0	0.0	-6916463.448159969;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4079026487214111	0.0	0.0	6916463.448159969;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.4079026487214111	0.0	0.0	-6916463.448159969;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.362169339054259	0.0	0.0	6916463.448159969;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.4215612663095952	0.0	0.0	6916463.448159969;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.4215612663095952	0.0	0.0	-6916463.448159969;	%DE1 76 PHS (pump mode)
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

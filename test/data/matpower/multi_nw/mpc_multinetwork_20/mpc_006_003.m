function mpc = mpc_006_003
%MPC_006_003	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 6 	Weight: 36
%	Time step: 3

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	-17.35472716536074	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	4231.021257198495	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	11	2	1201.1943858475786	0.0	0.0	0.0	11	1.0	0.0	230.0	1	1.1	0.9;
	12	2	-342.0921239367705	0.0	0.0	0.0	12	1.0	0.0	230.0	1	1.1	0.9;
	13	2	2166.7247027054254	0.0	0.0	0.0	13	1.0	0.0	230.0	1	1.1	0.9;
	14	2	5248.115005963613	0.0	0.0	0.0	14	1.0	0.0	230.0	1	1.1	0.9;
	15	2	231.89086285907902	0.0	0.0	0.0	15	1.0	0.0	230.0	1	1.1	0.9;
	16	2	2282.445122142996	0.0	0.0	0.0	16	1.0	0.0	230.0	1	1.1	0.9;
	17	2	1118.4047319127892	0.0	0.0	0.0	17	1.0	0.0	230.0	1	1.1	0.9;
	18	2	3090.8483548792656	0.0	0.0	0.0	18	1.0	0.0	230.0	1	1.1	0.9;
	19	2	1395.7957110662865	0.0	0.0	0.0	19	1.0	0.0	230.0	1	1.1	0.9;
	20	2	3869.123637064471	0.0	0.0	0.0	20	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4546.7522980816075	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2941.542823495351	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	3210.6488671302172	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	3879.1387364937277	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	4985.079874956537	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4061.2653095040732	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3038.0359201613837	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	3304.586695010812	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.5255	2.0436	2.8	8.0	8.0	469.0;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.750604815374423	0.0	0.0	1.0	1.0	276.0;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.8292451245865773	0.0	0.0	43.0	43.0	124.0;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	23.0	23.0	168.0;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.82	78.4	5.0	5.0	5.0	53.0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.10750000000002	78.4	2.8	2.0	2.0	53.0;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.716843749999999	2.0436	2.8	24.0	24.0	1172.0;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.1292490730129954	0.0	0.0	57.0	57.0	309.0;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	75.0	75.0	419.0;	%1 solar
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	174.25	78.4	5.0	2.0	2.0	74.0;	%10 CCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	149.3	78.4	2.8	1.0	1.0	74.0;	%10 OCGT
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.873419199999997	2.0436	2.8	50.0	50.0	1637.0;	%10 biomass
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.9155248407441183	0.0	0.0	20.0	20.0	963.0;	%10 offwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.4417809050615005	0.0	0.0	118.0	118.0	431.0;	%10 onwind
	11	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	86.0	86.0	585.0;	%10 solar
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	46.0;	%11 CCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.5	78.4	2.8	1.0	1.0	46.0;	%11 OCGT
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.822157142857145	2.0436	2.8	35.0	35.0	1023.0;	%11 biomass
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.98742248007162	0.0	0.0	7.0	7.0	602.0;	%11 offwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.126606318540895	0.0	0.0	86.0	86.0	270.0;	%11 onwind
	12	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	40.0	40.0	366.0;	%11 solar
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	58.0;	%12 OCGT
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.90428347826087	2.0436	2.8	69.0	69.0	1280.0;	%12 biomass
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.055407149480815	0.0	0.0	29.0	29.0	337.0;	%12 onwind
	13	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	166.0	166.0	457.0;	%12 solar
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.29199999999997	78.4	5.0	5.0	5.0	72.0;	%13 CCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	198.34283333333335	78.4	2.8	3.0	3.0	72.0;	%13 OCGT
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.778892954545455	2.0436	2.8	44.0	44.0	1597.0;	%13 biomass
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.148790199011488	0.0	0.0	45.0	45.0	421.0;	%13 onwind
	14	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	154.0	154.0	571.0;	%13 solar
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	5.0	1.0	1.0	28.0;	%14 CCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	28.0;	%14 OCGT
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.49785	2.0436	2.8	11.0	11.0	621.0;	%14 biomass
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.751599057472453	0.0	0.0	34.0	34.0	365.0;	%14 offwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.353496079888187	0.0	0.0	37.0	37.0	164.0;	%14 onwind
	15	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33.0	33.0	222.0;	%14 solar
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	142.0	78.4	5.0	2.0	2.0	75.0;	%15 CCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	75.0;	%15 OCGT
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.853806481481481	2.0436	2.8	54.0	54.0	1672.0;	%15 biomass
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.7821303105533657	0.0	0.0	46.0	46.0	984.0;	%15 offwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.508246538860008	0.0	0.0	137.0	137.0	440.0;	%15 onwind
	16	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56.0	56.0	597.0;	%15 solar
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	78.4	5.0	2.0	2.0	58.0;	%16 CCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	133.45	78.4	2.8	1.0	1.0	58.0;	%16 OCGT
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93894076923077	2.0436	2.8	52.0	52.0	1279.0;	%16 biomass
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.103898842297502	0.0	0.0	18.0	18.0	337.0;	%16 onwind
	17	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	173.0	173.0	457.0;	%16 solar
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	172.2	78.4	5.0	1.0	1.0	45.0;	%17 CCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	80.11	78.4	2.8	1.0	1.0	45.0;	%17 OCGT
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.823858333333332	2.0436	2.8	12.0	12.0	991.0;	%17 biomass
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.5418940087456585	0.0	0.0	72.0	72.0	261.0;	%17 onwind
	18	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	67.0	67.0	354.0;	%17 solar
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	215.0	78.4	5.0	2.0	2.0	44.0;	%18 CCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	44.0;	%18 OCGT
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.530666666666667	2.0436	2.8	18.0	18.0	974.0;	%18 biomass
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.882100965020366	0.0	0.0	53.0	53.0	257.0;	%18 onwind
	19	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	81.0	81.0	348.0;	%18 solar
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.9294117647059	78.4	5.0	17.0	17.0	70.0;	%19 CCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	168.36975	78.4	2.8	4.0	4.0	70.0;	%19 OCGT
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.898548538461538	2.0436	2.8	13.0	13.0	1562.0;	%19 biomass
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.057991255220416	0.0	0.0	49.0	49.0	411.0;	%19 onwind
	20	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	56.0	56.0	558.0;	%19 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	199.975	78.4	5.0	4.0	4.0	81.0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	138.3	78.4	2.8	1.0	1.0	81.0;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908959523809523	2.0436	2.8	63.0	63.0	1794.0;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.404160455000066	0.0	0.0	82.0	82.0	472.0;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	155.0	155.0	641.0;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	192.35999999999996	78.4	5.0	5.0	5.0	139.0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	78.4	78.4	2.8	1.0	1.0	139.0;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.850204846153847	2.0436	2.8	65.0	65.0	3109.0;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.337245658763658	0.0	0.0	350.0	350.0	1829.0;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.333857413876464	0.0	0.0	156.0	156.0	818.0;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	73.0	73.0	1111.0;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	206.15833333333333	78.4	5.0	12.0	12.0	58.0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	127.455	78.4	2.8	2.0	2.0	58.0;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.779668421052632	2.0436	2.8	38.0	38.0	1279.0;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8031995520920023	0.0	0.0	16.0	16.0	337.0;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	98.0	98.0	457.0;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2.0	2.0	35.0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	89.0	78.4	2.8	1.0	1.0	35.0;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.905931481481481	2.0436	2.8	27.0	27.0	779.0;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.548609146326785	0.0	0.0	17.0	17.0	205.0;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	90.0	90.0	278.0;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	212.0017272727273	78.4	5.0	11.0	11.0	73.0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	158.5	78.4	2.8	1.0	1.0	73.0;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.929266666666667	2.0436	2.8	45.0	45.0	1622.0;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.6198249852106286	0.0	0.0	70.0	70.0	427.0;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	91.0	91.0	580.0;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.0	78.4	5.0	2.0	2.0	76.0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	123.26000000000002	78.4	2.8	2.0	2.0	76.0;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.831483333333335	2.0436	2.8	54.0	54.0	1682.0;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.3167693295849823	0.0	0.0	109.0	109.0	443.0;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	111.0	111.0	601.0;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	184.222	78.4	5.0	5.0	5.0	101.0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	221.50666666666666	78.4	2.8	3.0	3.0	101.0;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.856493333333333	2.0436	2.8	60.0	60.0	2241.0;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.6488687661110744	0.0	0.0	113.0	113.0	590.0;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170.0	170.0	801.0;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	209.9742857142857	78.4	5.0	14.0	14.0	112.0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.4	78.4	2.8	1.0	1.0	112.0;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.998257142857144	2.0436	2.8	70.0	70.0	2491.0;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.8542786397604463	0.0	0.0	133.0	133.0	656.0;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	132.0	132.0	890.0;	%9 solar
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
	2	2741.232028703249	1096.4928114812997	3	0.0	536.5496239387843	0.0	0.0	10010990.352817608;	%0 biomass
	2	0.0	0.0	3	0.0	0.888224483995433	0.0	0.0	6820180.090177206;	%0 offwind
	2	0.0	0.0	3	0.0	0.8675093602255116	0.0	0.0	3459091.9687187863;	%0 onwind
	2	0.0	0.0	3	0.0	0.6861189854382682	0.0	0.0	1281674.5647769603;	%0 solar
	2	319784.5005205832	127913.80020823325	3	0.0	1631.5535740846076	0.0	0.0	3564998.254553488;	%1 CCGT
	2	390347.9161637284	156139.16646549135	3	0.0	1991.571000835349	0.0	0.0	1623435.8945444098;	%1 OCGT
	2	2741.122521953878	1096.4490087815511	3	0.0	536.5281898520019	0.0	0.0	10010990.352817606;	%1 biomass
	2	0.0	0.0	3	0.0	0.8848805037083349	0.0	0.0	3459091.9687187863;	%1 onwind
	2	0.0	0.0	3	0.0	0.7178316710257875	0.0	0.0	1281674.5647769603;	%1 solar
	2	280197.9775864519	112079.19103458076	3	0.0	1429.5815182982237	0.0	0.0	3564998.254553488;	%10 CCGT
	2	406447.4677879188	162578.98711516752	3	0.0	2073.7115703465242	0.0	0.0	1623435.8945444098;	%10 OCGT
	2	2741.1729475682196	1096.4691790272877	3	0.0	536.5380598097904	0.0	0.0	10010990.352817606;	%10 biomass
	2	0.0	0.0	3	0.0	0.8889511874905132	0.0	0.0	6733080.78488266;	%10 offwind
	2	0.0	0.0	3	0.0	0.8927895337030355	0.0	0.0	3459091.968718787;	%10 onwind
	2	0.0	0.0	3	0.0	0.7136020083283926	0.0	0.0	1281674.5647769603;	%10 solar
	2	261452.88451692156	104581.15380676862	3	0.0	1333.9432883516406	0.0	0.0	3564998.254553488;	%11 CCGT
	2	388126.9757205369	155250.79028821475	3	0.0	1980.2396720435554	0.0	0.0	1623435.8945444098;	%11 OCGT
	2	2741.2335711601677	1096.4934284640672	3	0.0	536.5499258485354	0.0	0.0	10010990.352817604;	%11 biomass
	2	0.0	0.0	3	0.0	0.9113119523434836	0.0	0.0	6694538.863358395;	%11 offwind
	2	0.0	0.0	3	0.0	0.8678451707189306	0.0	0.0	3459091.9687187867;	%11 onwind
	2	0.0	0.0	3	0.0	0.7234344921369239	0.0	0.0	1281674.5647769603;	%11 solar
	2	386035.33073458896	154414.13229383557	3	0.0	1969.5680139519843	0.0	0.0	1623435.8945444098;	%12 OCGT
	2	2741.181472237433	1096.4725888949733	3	0.0	536.5397283690415	0.0	0.0	10010990.352817606;	%12 biomass
	2	0.0	0.0	3	0.0	0.9005743048745756	0.0	0.0	3459091.9687187867;	%12 onwind
	2	0.0	0.0	3	0.0	0.7203638607337804	0.0	0.0	1281674.5647769603;	%12 solar
	2	290226.6315195534	116090.65260782136	3	0.0	1480.7481199977212	0.0	0.0	3564998.254553488;	%13 CCGT
	2	403718.0018988323	161487.20075953292	3	0.0	2059.785723973634	0.0	0.0	1623435.8945444098;	%13 OCGT
	2	2741.1563098621436	1096.4625239448574	3	0.0	536.5348032613317	0.0	0.0	10010990.352817606;	%13 biomass
	2	0.0	0.0	3	0.0	0.9154525948192558	0.0	0.0	3459091.9687187867;	%13 onwind
	2	0.0	0.0	3	0.0	0.7287711180840533	0.0	0.0	1281674.5647769603;	%13 solar
	2	299077.0560698286	119630.82242793146	3	0.0	1525.9033472950441	0.0	0.0	3564998.254553488;	%14 CCGT
	2	403428.9647241115	161371.5858896446	3	0.0	2058.311044510773	0.0	0.0	1623435.8945444098;	%14 OCGT
	2	2741.0439537927946	1096.417581517118	3	0.0	536.5128114685447	0.0	0.0	10010990.352817606;	%14 biomass
	2	0.0	0.0	3	0.0	0.8900895727961411	0.0	0.0	7050536.01427284;	%14 offwind
	2	0.0	0.0	3	0.0	0.8812800900780561	0.0	0.0	3459091.968718787;	%14 onwind
	2	0.0	0.0	3	0.0	0.7045151104027109	0.0	0.0	1281674.5647769603;	%14 solar
	2	268597.5661364501	107439.02645458003	3	0.0	1370.395745594133	0.0	0.0	3564998.254553488;	%15 CCGT
	2	376733.01669747365	150693.20667898946	3	0.0	1922.107228048335	0.0	0.0	1623435.8945444098;	%15 OCGT
	2	2741.187558679373	1096.4750234717492	3	0.0	536.5409196867045	0.0	0.0	10010990.352817606;	%15 biomass
	2	0.0	0.0	3	0.0	0.9011508992295109	0.0	0.0	7357937.603524846;	%15 offwind
	2	0.0	0.0	3	0.0	0.8920342906062702	0.0	0.0	3459091.968718787;	%15 onwind
	2	0.0	0.0	3	0.0	0.6962192497405517	0.0	0.0	1281674.5647769603;	%15 solar
	2	281161.17268683005	112464.46907473201	3	0.0	1434.4957790144388	0.0	0.0	3564998.254553488;	%16 CCGT
	2	389020.9079580178	155608.36318320708	3	0.0	1984.8005508062129	0.0	0.0	1623435.89454441;	%16 OCGT
	2	2741.121505043071	1096.4486020172285	3	0.0	536.5279908089784	0.0	0.0	10010990.352817608;	%16 biomass
	2	0.0	0.0	3	0.0	0.890523474937217	0.0	0.0	3459091.9687187867;	%16 onwind
	2	0.0	0.0	3	0.0	0.7161234305804126	0.0	0.0	1281674.5647769605;	%16 solar
	2	289070.3361621399	115628.13446485596	3	0.0	1474.8486538884686	0.0	0.0	3564998.254553488;	%17 CCGT
	2	379232.9592486978	151693.18369947912	3	0.0	1934.8620369831517	0.0	0.0	1623435.8945444098;	%17 OCGT
	2	2741.152459839751	1096.4609839359005	3	0.0	536.5340496848211	0.0	0.0	10010990.352817606;	%17 biomass
	2	0.0	0.0	3	0.0	0.9177161030789063	0.0	0.0	3459091.968718787;	%17 onwind
	2	0.0	0.0	3	0.0	0.7338236927614148	0.0	0.0	1281674.5647769605;	%17 solar
	2	303010.8616640877	121204.34466563509	3	0.0	1545.9737840004473	0.0	0.0	3564998.254553488;	%18 CCGT
	2	382184.3803223054	152873.75212892218	3	0.0	1949.9203077668644	0.0	0.0	1623435.8945444098;	%18 OCGT
	2	2741.20207444201	1096.480829776804	3	0.0	536.5437609007653	0.0	0.0	10010990.352817606;	%18 biomass
	2	0.0	0.0	3	0.0	0.8967454302943216	0.0	0.0	3459091.9687187867;	%18 onwind
	2	0.0	0.0	3	0.0	0.7199796322909876	0.0	0.0	1281674.5647769603;	%18 solar
	2	283925.40429754247	113570.16171901699	3	0.0	1448.5990015180737	0.0	0.0	3564998.254553488;	%19 CCGT
	2	397805.1526976923	159122.0610790769	3	0.0	2029.6181260086337	0.0	0.0	1623435.8945444098;	%19 OCGT
	2	2741.13891181217	1096.455564724868	3	0.0	536.5313978884654	0.0	0.0	10010990.352817606;	%19 biomass
	2	0.0	0.0	3	0.0	0.9149052910129215	0.0	0.0	3459091.9687187863;	%19 onwind
	2	0.0	0.0	3	0.0	0.7156498977125215	0.0	0.0	1281674.5647769603;	%19 solar
	2	303968.9495466029	121587.57981864116	3	0.0	1550.8619874826677	0.0	0.0	3564998.254553488;	%2 CCGT
	2	401732.8730848537	160693.14923394148	3	0.0	2049.6575157390494	0.0	0.0	1623435.89454441;	%2 OCGT
	2	2741.1255060327603	1096.450202413104	3	0.0	536.528773934774	0.0	0.0	10010990.352817606;	%2 biomass
	2	0.0	0.0	3	0.0	0.9016331201408072	0.0	0.0	3459091.9687187867;	%2 onwind
	2	0.0	0.0	3	0.0	0.7179486594564432	0.0	0.0	1281674.5647769603;	%2 solar
	2	331485.17599425855	132594.07039770344	3	0.0	1691.2508979298907	0.0	0.0	3564998.254553488;	%3 CCGT
	2	412342.841007491	164937.1364029964	3	0.0	2103.7900051402603	0.0	0.0	1623435.8945444098;	%3 OCGT
	2	2741.0454328974533	1096.4181731589813	3	0.0	536.5131009781667	0.0	0.0	10010990.352817606;	%3 biomass
	2	0.0	0.0	3	0.0	0.8918750555367175	0.0	0.0	7325177.817049337;	%3 offwind
	2	0.0	0.0	3	0.0	0.9175290142773271	0.0	0.0	3459091.968718787;	%3 onwind
	2	0.0	0.0	3	0.0	0.7189254862017672	0.0	0.0	1281674.5647769603;	%3 solar
	2	288959.9509182839	115583.98036731356	3	0.0	1474.2854638687954	0.0	0.0	3564998.254553488;	%4 CCGT
	2	395177.39466730325	158070.95786692132	3	0.0	2016.2111972821594	0.0	0.0	1623435.8945444098;	%4 OCGT
	2	2741.208240462323	1096.483296184929	3	0.0	536.5449677945435	0.0	0.0	10010990.352817606;	%4 biomass
	2	0.0	0.0	3	0.0	0.8974730745617018	0.0	0.0	3459091.9687187867;	%4 onwind
	2	0.0	0.0	3	0.0	0.7036033915681454	0.0	0.0	1281674.5647769603;	%4 solar
	2	296799.5295453994	118719.81181815975	3	0.0	1514.2833140071398	0.0	0.0	3564998.254553488;	%5 CCGT
	2	388056.7171060436	155222.68684241743	3	0.0	1979.8812097247119	0.0	0.0	1623435.8945444098;	%5 OCGT
	2	2741.1717690091778	1096.4687076036712	3	0.0	536.5378291268698	0.0	0.0	10010990.352817604;	%5 biomass
	2	0.0	0.0	3	0.0	0.8865621856068566	0.0	0.0	3459091.9687187867;	%5 onwind
	2	0.0	0.0	3	0.0	0.7206121964851012	0.0	0.0	1281674.5647769603;	%5 solar
	2	291009.0254825078	116403.61019300311	3	0.0	1484.7399259311621	0.0	0.0	3564998.254553488;	%6 CCGT
	2	401917.0255072881	160766.81020291525	3	0.0	2050.5970689147352	0.0	0.0	1623435.89454441;	%6 OCGT
	2	2741.217223163506	1096.4868892654024	3	0.0	536.5467260057753	0.0	0.0	10010990.352817606;	%6 biomass
	2	0.0	0.0	3	0.0	0.9003488870087848	0.0	0.0	3459091.9687187867;	%6 onwind
	2	0.0	0.0	3	0.0	0.7350788485422832	0.0	0.0	1281674.5647769603;	%6 solar
	2	275834.4258555436	110333.77034221742	3	0.0	1407.3184992629774	0.0	0.0	3564998.254553488;	%7 CCGT
	2	400198.7463575459	160079.49854301836	3	0.0	2041.8303385589074	0.0	0.0	1623435.8945444098;	%7 OCGT
	2	2741.21428769085	1096.48571507634	3	0.0	536.5461514368468	0.0	0.0	10010990.352817604;	%7 biomass
	2	0.0	0.0	3	0.0	0.9101968393491763	0.0	0.0	3459091.9687187867;	%7 onwind
	2	0.0	0.0	3	0.0	0.7182037471582787	0.0	0.0	1281674.5647769605;	%7 solar
	2	311869.8034196855	124747.9213678742	3	0.0	1591.1724664269666	0.0	0.0	3564998.254553488;	%8 CCGT
	2	408979.9537951865	163591.9815180746	3	0.0	2086.63241732238	0.0	0.0	1623435.8945444098;	%8 OCGT
	2	2741.2194985647106	1096.487799425884	3	0.0	536.5471713769251	0.0	0.0	10010990.352817606;	%8 biomass
	2	0.0	0.0	3	0.0	0.9200762907491696	0.0	0.0	3459091.9687187867;	%8 onwind
	2	0.0	0.0	3	0.0	0.7172869998323175	0.0	0.0	1281674.5647769603;	%8 solar
	2	299056.70019625616	119622.68007850247	3	0.0	1525.7994907972252	0.0	0.0	3564998.254553488;	%9 CCGT
	2	403509.83077638416	161403.93231055368	3	0.0	2058.723626410123	0.0	0.0	1623435.8945444098;	%9 OCGT
	2	2741.111918311479	1096.4447673245916	3	0.0	536.5261143690504	0.0	0.0	10010990.352817606;	%9 biomass
	2	0.0	0.0	3	0.0	0.8948063584329272	0.0	0.0	3459091.9687187863;	%9 onwind
	2	0.0	0.0	3	0.0	0.725646810536715	0.0	0.0	1281674.5647769603;	%9 solar
	2	0.0	0.0	3	0.0	0.34805112935570215	0.0	0.0	6384427.7983015105;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.34805112935570215	0.0	0.0	-6384427.7983015105;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3791950700194782	0.0	0.0	6384427.7983015105;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.3791950700194782	0.0	0.0	-6384427.7983015105;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39566972183147336	0.0	0.0	6384427.7983015105;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.39566972183147336	0.0	0.0	-6384427.7983015105;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3304239874391053	0.0	0.0	6384427.7983015105;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.3304239874391053	0.0	0.0	-6384427.7983015105;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39310872116804446	0.0	0.0	6384427.7983015105;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.39310872116804446	0.0	0.0	-6384427.7983015105;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.378170855844992	0.0	0.0	6384427.7983015105;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.378170855844992	0.0	0.0	-6384427.7983015105;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3599684561766089	0.0	0.0	6384427.7983015105;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.3599684561766089	0.0	0.0	-6384427.7983015105;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.36539495935741817	0.0	0.0	6384427.7983015105;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.361282317054359	0.0	0.0	6384427.7983015105;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.361282317054359	0.0	0.0	-6384427.7983015105;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3518076716409871	0.0	0.0	6384427.7983015105;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.3518076716409871	0.0	0.0	-6384427.7983015105;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3453413335152406	0.0	0.0	6384427.7983015105;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.3453413335152406	0.0	0.0	-6384427.7983015105;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3307149104559095	0.0	0.0	6384427.7983015105;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.3307149104559095	0.0	0.0	-6384427.7983015105;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3647682750306956	0.0	0.0	6384427.7983015105;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.3647682750306956	0.0	0.0	-6384427.7983015105;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33483945325334463	0.0	0.0	6384427.7983015105;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.33483945325334463	0.0	0.0	-6384427.7983015105;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3368539160946931	0.0	0.0	6384427.7983015105;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.3368539160946931	0.0	0.0	-6384427.7983015105;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3903708631833174	0.0	0.0	6384427.7983015105;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.3903708631833174	0.0	0.0	-6384427.7983015105;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.35735422688248913	0.0	0.0	6384427.7983015105;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.35735422688248913	0.0	0.0	-6384427.7983015105;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3714009329509926	0.0	0.0	6384427.7983015105;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.3714009329509926	0.0	0.0	-6384427.7983015105;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37652552189668714	0.0	0.0	6384427.7983015105;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.37652552189668714	0.0	0.0	-6384427.7983015105;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33431015912700834	0.0	0.0	6384427.7983015105;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.38913347659347247	0.0	0.0	6384427.7983015105;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.38913347659347247	0.0	0.0	-6384427.7983015105;	%DE1 76 PHS (pump mode)
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

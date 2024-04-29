function mpc = mpc_012_000
%MPC_012_000	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 12 	Weight: 37
%	Time step: 0

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 37;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	838.3872698651531	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	7769.305893100084	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	3639.6361842784936	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	1806.087509492398	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	6510.433436432588	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	8969.360527943423	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	9112.42695998424	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4546.016456965471	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	2840.017194346498	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	7003.0022537809045	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.569940114939909	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.96346696892241	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8796711239152761	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8918916640796716	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.9627743859062043	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.363438388535558	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8194361182840926	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.2629946575723028	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.5435430332008933	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.8584341696303468	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.7365810826843024	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.490820942499961	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.860136599542438	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	221	221	2191	28.0	5;	%9 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1776.8	266.52	3.5	1	1	1	1776.8	7;	%DE1 0 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-266.52	-1776.8	3.5	1	1	1	-1776.8	8;	%DE1 0 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	1294.0	194.1	3.5	1	1	1	1294.0	7;	%DE1 12 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-194.1	-1294.0	3.5	1	1	1	-1294.0	8;	%DE1 12 PHS (pump mode)
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
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	167.1	25.064999999999998	3.5	1	1	1	167.1	7;	%DE1 27 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-25.064999999999998	-167.1	3.5	1	1	1	-167.1	8;	%DE1 27 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	109.0	16.349999999999998	3.5	1	1	1	109.0	7;	%DE1 31 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-16.349999999999998	-109.0	3.5	1	1	1	-109.0	8;	%DE1 31 PHS (pump mode)
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	480.0	72.0	3.5	1	1	1	480.0	7;	%DE1 33 PHS
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	-72.0	-480.0	3.5	1	1	1	-480.0	8;	%DE1 33 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	92.0	13.799999999999999	3.5	1	1	1	92.0	7;	%DE1 36 PHS
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.799999999999999	-92.0	3.5	1	1	1	-92.0	8;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	278.0	41.699999999999996	3.5	1	1	1	278.0	7;	%DE1 43 PHS
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	-41.699999999999996	-278.0	3.5	1	1	1	-278.0	8;	%DE1 43 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	164.0	24.599999999999998	3.5	1	1	1	164.0	7;	%DE1 49 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.599999999999998	-164.0	3.5	1	1	1	-164.0	8;	%DE1 49 PHS (pump mode)
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	79.7	11.955	3.5	1	1	1	79.7	7;	%DE1 50 PHS
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	-11.955	-79.7	3.5	1	1	1	-79.7	8;	%DE1 50 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	399.8	59.97	3.5	1	1	1	399.8	7;	%DE1 52 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-59.97	-399.8	3.5	1	1	1	-399.8	8;	%DE1 52 PHS (pump mode)
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	162.0	24.3	3.5	1	1	1	162.0	7;	%DE1 54 PHS
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	-24.3	-162.0	3.5	1	1	1	-162.0	8;	%DE1 54 PHS (pump mode)
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	90.0	13.5	3.5	1	1	1	90.0	7;	%DE1 61 PHS
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	-13.5	-90.0	3.5	1	1	1	-90.0	8;	%DE1 61 PHS (pump mode)
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	1060.0	159.0	3.5	1	1	1	1060.0	7;	%DE1 7 PHS
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	-159.0	-1060.0	3.5	1	1	1	-1060.0	8;	%DE1 7 PHS (pump mode)
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	45.5	6.825	3.5	1	1	1	45.5	6;	%DE1 71 hydro
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.0	18.0	3.5	1	1	1	120.0	7;	%DE1 76 PHS
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	-18.0	-120.0	3.5	1	1	1	-120.0	8;	%DE1 76 PHS (pump mode)
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0	investment
mpc.gencost = [
	2	0	0	3	0.0	1574.3801622789827	0.0	0.0	93691.98862804195;	%0 CCGT
	2	0	0	3	0.0	2123.8659362864732	0.0	0.0	42665.64147001544;	%0 OCGT
	2	0	0	3	0.0	551.4434953653879	0.0	0.0	11745.517600654852;	%0 biomass
	2	0	0	3	0.0	0.9142061529998328	0.0	0.0	13816.073021315466;	%0 offwind
	2	0	0	3	0.0	0.9093578359882806	0.0	0.0	15422.00440137815;	%0 onwind
	2	0	0	3	0.0	0.7276997880946361	0.0	0.0	4210.473265819795;	%0 solar
	2	0	0	3	0.0	1581.4524739440185	0.0	0.0	93691.98862804195;	%1 CCGT
	2	0	0	3	0.0	2034.0658700193053	0.0	0.0	42665.64147001544;	%1 OCGT
	2	0	0	3	0.0	551.4318830603366	0.0	0.0	11745.517600654852;	%1 biomass
	2	0	0	3	0.0	0.9255393067802922	0.0	0.0	15422.00440137815;	%1 onwind
	2	0	0	3	0.0	0.7448340626906325	0.0	0.0	4210.473265819795;	%1 solar
	2	0	0	3	0.0	1507.6029534726197	0.0	0.0	93691.98862804195;	%2 CCGT
	2	0	0	3	0.0	2089.979016642353	0.0	0.0	42665.64147001544;	%2 OCGT
	2	0	0	3	0.0	551.4273360074181	0.0	0.0	11745.517600654852;	%2 biomass
	2	0	0	3	0.0	0.9291815675300901	0.0	0.0	15422.00440137815;	%2 onwind
	2	0	0	3	0.0	0.7385495765652935	0.0	0.0	4210.473265819795;	%2 solar
	2	0	0	3	0.0	1395.9739522497828	0.0	0.0	93691.98862804195;	%3 CCGT
	2	0	0	3	0.0	2005.3727125472212	0.0	0.0	42665.64147001544;	%3 OCGT
	2	0	0	3	0.0	551.4430352934683	0.0	0.0	11745.517600654852;	%3 biomass
	2	0	0	3	0.0	0.9303601904882973	0.0	0.0	14146.47269657296;	%3 offwind
	2	0	0	3	0.0	0.9137113974759555	0.0	0.0	15422.004401378152;	%3 onwind
	2	0	0	3	0.0	0.7344617253122686	0.0	0.0	4210.473265819795;	%3 solar
	2	0	0	3	0.0	1494.7903609261064	0.0	0.0	93691.98862804195;	%4 CCGT
	2	0	0	3	0.0	2063.040029159337	0.0	0.0	42665.64147001544;	%4 OCGT
	2	0	0	3	0.0	551.4436329245079	0.0	0.0	11745.517600654852;	%4 biomass
	2	0	0	3	0.0	0.9161491646480641	0.0	0.0	15422.00440137815;	%4 onwind
	2	0	0	3	0.0	0.7235498384491076	0.0	0.0	4210.473265819795;	%4 solar
	2	0	0	3	0.0	1556.346739396227	0.0	0.0	93691.98862804195;	%5 CCGT
	2	0	0	3	0.0	2072.6051926071527	0.0	0.0	42665.64147001544;	%5 OCGT
	2	0	0	3	0.0	551.4462877835542	0.0	0.0	11745.517600654852;	%5 biomass
	2	0	0	3	0.0	0.9256517348450279	0.0	0.0	15422.00440137815;	%5 onwind
	2	0	0	3	0.0	0.7489381176225196	0.0	0.0	4210.473265819795;	%5 solar
	2	0	0	3	0.0	1480.2976490118815	0.0	0.0	93691.98862804195;	%6 CCGT
	2	0	0	3	0.0	2090.963534492936	0.0	0.0	42665.64147001544;	%6 OCGT
	2	0	0	3	0.0	551.438565342188	0.0	0.0	11745.517600654852;	%6 biomass
	2	0	0	3	0.0	0.9386322506284962	0.0	0.0	15422.00440137815;	%6 onwind
	2	0	0	3	0.0	0.7474630936006486	0.0	0.0	4210.473265819795;	%6 solar
	2	0	0	3	0.0	1504.38834323546	0.0	0.0	93691.98862804195;	%7 CCGT
	2	0	0	3	0.0	2104.9150419257426	0.0	0.0	42665.64147001544;	%7 OCGT
	2	0	0	3	0.0	551.4457928413836	0.0	0.0	11745.517600654852;	%7 biomass
	2	0	0	3	0.0	0.9313707311191337	0.0	0.0	15422.00440137815;	%7 onwind
	2	0	0	3	0.0	0.734317110197582	0.0	0.0	4210.473265819795;	%7 solar
	2	0	0	3	0.0	1617.951376086939	0.0	0.0	93691.98862804195;	%8 CCGT
	2	0	0	3	0.0	2116.492495283812	0.0	0.0	42665.64147001544;	%8 OCGT
	2	0	0	3	0.0	551.450091066859	0.0	0.0	11745.517600654852;	%8 biomass
	2	0	0	3	0.0	0.9376409855215597	0.0	0.0	15422.00440137815;	%8 onwind
	2	0	0	3	0.0	0.7381341146514633	0.0	0.0	4210.473265819795;	%8 solar
	2	0	0	3	0.0	1629.3922893550496	0.0	0.0	93691.98862804195;	%9 CCGT
	2	0	0	3	0.0	2124.4581516075323	0.0	0.0	42665.64147001544;	%9 OCGT
	2	0	0	3	0.0	551.4283886130971	0.0	0.0	11745.517600654852;	%9 biomass
	2	0	0	3	0.0	0.9166493626349597	0.0	0.0	14610.403497349165;	%9 offwind
	2	0	0	3	0.0	0.9241983826819631	0.0	0.0	15422.00440137815;	%9 onwind
	2	0	0	3	0.0	0.7396813984429169	0.0	0.0	4210.473265819795;	%9 solar
	2	0.0	0.0	3	0.0	0.35771921628224945	0.0	0.0	749.0608464521686;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.2682894122116871	0.0	0.0	749.0608464521686;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3897282664089082	0.0	0.0	749.0608464521686;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.29229619980668115	0.0	0.0	749.0608464521686;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40666054743790314	0.0	0.0	749.0608464521686;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.3049954105784274	0.0	0.0	749.0608464521686;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.339602431534636	0.0	0.0	749.0608464521686;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.254701823650977	0.0	0.0	749.0608464521686;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.40402840786715677	0.0	0.0	749.0608464521686;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.3030213059003676	0.0	0.0	749.0608464521686;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3886756018406862	0.0	0.0	749.0608464521686;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.29150670138051465	0.0	0.0	749.0608464521686;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.36996757995929247	0.0	0.0	749.0608464521686;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.27747568496946934	0.0	0.0	749.0608464521686;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3755448193395686	0.0	0.0	749.0608464521686;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.3713179369725356	0.0	0.0	749.0608464521686;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.2784884527294017	0.0	0.0	749.0608464521686;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.36158010696434784	0.0	0.0	749.0608464521686;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.2711850802232609	0.0	0.0	749.0608464521686;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3549341483351084	0.0	0.0	749.0608464521686;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.2662006112513313	0.0	0.0	749.0608464521686;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3399014357463514	0.0	0.0	749.0608464521686;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.25492607680976354	0.0	0.0	749.0608464521686;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37490072711488154	0.0	0.0	749.0608464521686;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.2811755453361612	0.0	0.0	749.0608464521686;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34414054917704867	0.0	0.0	749.0608464521686;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.2581054118827865	0.0	0.0	749.0608464521686;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3462109693195457	0.0	0.0	749.0608464521686;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.2596582269896593	0.0	0.0	749.0608464521686;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4012144982717429	0.0	0.0	749.0608464521686;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.30091087370380715	0.0	0.0	749.0608464521686;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3672807331847805	0.0	0.0	749.0608464521686;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.27546054988858537	0.0	0.0	749.0608464521686;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3817176255329646	0.0	0.0	749.0608464521686;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.28628821914972347	0.0	0.0	749.0608464521686;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.38698456417159516	0.0	0.0	749.0608464521686;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.29023842312869635	0.0	0.0	749.0608464521686;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3435965524360919	0.0	0.0	749.0608464521686;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.39994273983218004	0.0	0.0	749.0608464521686;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.29995705487413504	0.0	0.0	749.0608464521686;	%DE1 76 PHS (pump mode)
];
%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0.002346193481236678	0.019238786546140762	1.958387880013024e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	1	8	0.0023321905176023603	0.019123962244339356	3.236859740616215e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	1	9	0.0019445334507215365	0.0159451742959166	2.6988284162105606e-05	4379.317261857151	4379.317261857151	4379.317261857151	0	0	1	-60	60;
	2	3	0.002476449448812836	0.020306885480265262	2.067113656570062e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	2	6	0.00032235809806776487	0.002643336404155672	9.943293073955649e-05	20645.352805897994	20645.352805897994	20645.352805897994	0	0	1	-60	60;
	2	7	0.00018528201517929854	0.001519312524470248	8.729978747518073e-05	25516.226086943192	25516.226086943192	25516.226086943192	0	0	1	-60	60;
	2	8	0.0027328436533095644	0.02240931795713843	2.2811281044864822e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	2	10	0.000789978277571624	0.006477821876087317	8.077670377089266e-05	11886.718282183694	11886.718282183694	11886.718282183694	0	0	1	-60	60;
	3	5	0.0007685374454869242	0.006302007052992778	4.9916538212103785e-05	9473.625097078733	9473.625097078733	9473.625097078733	0	0	1	-60	60;
	3	6	0.0011358075415830791	0.00931362184098125	4.969480021174222e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	3	8	0.0007808106224918078	0.006402647104432824	3.4162678506126354e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	3	9	0.0009715277760101635	0.007966527763283342	4.250709469818385e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
	4	8	0.002553846799021281	0.0209415437519745	2.1317178905369238e-05	3396.205223481056	3396.205223481056	3396.205223481056	0	0	1	-60	60;
	4	10	0.0007315818716301728	0.0059989713473674175	4.061465629624203e-05	8758.634523714301	8758.634523714301	8758.634523714301	0	0	1	-60	60;
	5	6	0.0007449955480032163	0.0061089634936263735	4.838749101703001e-05	9473.625097078733	9473.625097078733	9473.625097078733	0	0	1	-60	60;
	7	10	0.00047207592830423385	0.0038710226120947182	6.304733408087732e-05	13584.820893924221	13584.820893924221	13584.820893924221	0	0	1	-60	60;
	8	9	0.0005598156111912242	0.00459048801176804	5.056298758014288e-05	11171.727708819262	11171.727708819262	11171.727708819262	0	0	1	-60	60;
	8	10	0.0008717140618085548	0.0071480553068301485	3.813996170774132e-05	7775.522485338206	7775.522485338206	7775.522485338206	0	0	1	-60	60;
];

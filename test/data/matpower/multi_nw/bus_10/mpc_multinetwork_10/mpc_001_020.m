function mpc = mpc_001_020
%MPC_001_020	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 1 	Weight: 30
%	Time step: 20

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 30;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	943.0910162502106	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	8434.050252267594	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4594.399385225968	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2790.163381363855	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	7151.331840569103	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	9059.623230141418	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	8843.967887631994	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4738.879545073967	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	3623.5075711736786	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	8097.564922824432	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.76303977756134	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	15.122057872539756	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.137068043281147	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.7346860209411545	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.193071445079315	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	14.47940434532662	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.834352467841429	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.544995645151684	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.924251152072244	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.098466895591107	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.634861532199414	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.505590934006213	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.292850190100086	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	250198.7933567681	100079.51734270723	3	0.0	1276.524455901878	0.0	0.0	339.1360592231248;	%0 CCGT
	2	337522.47852336394	135008.99140934556	3	0.0	1722.0534618538973	0.0	0.0	154.43644354494006;	%0 OCGT
	2	2284.3174198554866	913.7269679421947	3	0.0	447.11634759355775	0.0	0.0	952.339264917961;	%0 biomass
	2	0.0	0.0	3	0.0	0.7412482321620265	0.0	0.0	658.9541981549506;	%0 offwind
	2	0.0	0.0	3	0.0	0.7373171643148221	0.0	0.0	329.06126034235035;	%0 onwind
	2	0.0	0.0	3	0.0	0.590026855211867	0.0	0.0	121.92490151987825;	%0 solar
	2	251322.7174808332	100529.0869923333	3	0.0	1282.2587626573122	0.0	0.0	339.1360592231248;	%1 CCGT
	2	323251.5490733383	129300.61962933531	3	0.0	1649.2425973129502	0.0	0.0	154.43644354494006;	%1 OCGT
	2	2284.269316666427	913.7077266665708	3	0.0	447.10693221108374	0.0	0.0	952.339264917961;	%1 biomass
	2	0.0	0.0	3	0.0	0.7504372757678045	0.0	0.0	329.06126034235035;	%1 onwind
	2	0.0	0.0	3	0.0	0.603919510289702	0.0	0.0	121.92490151987825;	%1 solar
	2	239586.63152483798	95834.65260993519	3	0.0	1222.380773085908	0.0	0.0	339.1360592231248;	%2 CCGT
	2	332137.205888028	132854.88235521122	3	0.0	1694.5775810613673	0.0	0.0	154.43644354494006;	%2 OCGT
	2	2284.2504808069452	913.7001923227782	3	0.0	447.1032454114201	0.0	0.0	952.3392649179608;	%2 biomass
	2	0.0	0.0	3	0.0	0.7533904601595325	0.0	0.0	329.06126034235035;	%2 onwind
	2	0.0	0.0	3	0.0	0.5988239809988867	0.0	0.0	121.92490151987825;	%2 solar
	2	221846.67133050604	88738.66853220241	3	0.0	1131.8707720944185	0.0	0.0	339.1360592231248;	%3 CCGT
	2	318691.6635075044	127476.66540300177	3	0.0	1625.9778750382875	0.0	0.0	154.43644354494006;	%3 OCGT
	2	2284.3155140386452	913.7262056154581	3	0.0	447.1159745622715	0.0	0.0	952.339264917961;	%3 biomass
	2	0.0	0.0	3	0.0	0.7543461003959168	0.0	0.0	674.712529248313;	%3 offwind
	2	0.0	0.0	3	0.0	0.7408470790345586	0.0	0.0	329.0612603423504;	%3 onwind
	2	0.0	0.0	3	0.0	0.5955095070099475	0.0	0.0	121.92490151987825;	%3 solar
	2	237550.46816879747	95020.18726751898	3	0.0	1211.9921845346807	0.0	0.0	339.1360592231248;	%4 CCGT
	2	327856.0911204568	131142.4364481827	3	0.0	1672.7351587778408	0.0	0.0	154.43644354494006;	%4 OCGT
	2	2284.3179896848474	913.7271958739389	3	0.0	447.11645912797945	0.0	0.0	952.339264917961;	%4 biomass
	2	0.0	0.0	3	0.0	0.7428236470119438	0.0	0.0	329.06126034235035;	%4 onwind
	2	0.0	0.0	3	0.0	0.5866620311749521	0.0	0.0	121.92490151987825;	%4 solar
	2	247332.94128783283	98933.17651513313	3	0.0	1261.9027616726164	0.0	0.0	339.1360592231248;	%5 CCGT
	2	329376.17655486637	131750.47062194656	3	0.0	1680.4906967085021	0.0	0.0	154.43644354494006;	%5 OCGT
	2	2284.3289872590635	913.7315949036255	3	0.0	447.11861171639526	0.0	0.0	952.339264917961;	%5 biomass
	2	0.0	0.0	3	0.0	0.7505284336581308	0.0	0.0	329.06126034235035;	%5 onwind
	2	0.0	0.0	3	0.0	0.6072471223966375	0.0	0.0	121.92490151987825;	%5 solar
	2	235247.30205918552	94098.9208236742	3	0.0	1200.2413370366608	0.0	0.0	339.1360592231248;	%6 CCGT
	2	332293.6644004991	132917.46576019964	3	0.0	1695.3758387780563	0.0	0.0	154.43644354494006;	%6 OCGT
	2	2284.296997567491	913.7187990269963	3	0.0	447.11235027744976	0.0	0.0	952.339264917961;	%6 biomass
	2	0.0	0.0	3	0.0	0.7610531761852671	0.0	0.0	329.06126034235035;	%6 onwind
	2	0.0	0.0	3	0.0	0.6060511569734989	0.0	0.0	121.92490151987822;	%6 solar
	2	239075.76914120285	95630.30765648115	3	0.0	1219.7743323530758	0.0	0.0	339.1360592231248;	%7 CCGT
	2	334510.8228790099	133804.32915160398	3	0.0	1706.6878718316832	0.0	0.0	154.43644354494006;	%7 OCGT
	2	2284.326936994564	913.7307747978256	3	0.0	447.11821041193264	0.0	0.0	952.339264917961;	%7 biomass
	2	0.0	0.0	3	0.0	0.7551654576641624	0.0	0.0	329.06126034235035;	%7 onwind
	2	0.0	0.0	3	0.0	0.595392251511553	0.0	0.0	121.92490151987825;	%7 solar
	2	257123.08355111355	102849.23342044544	3	0.0	1311.8524670975182	0.0	0.0	339.1360592231248;	%8 CCGT
	2	336350.69925050857	134540.27970020342	3	0.0	1716.074996176064	0.0	0.0	154.43644354494006;	%8 OCGT
	2	2284.3447421031756	913.7378968412702	3	0.0	447.12169545961547	0.0	0.0	952.339264917961;	%8 biomass
	2	0.0	0.0	3	0.0	0.7602494477201835	0.0	0.0	329.06126034235035;	%8 onwind
	2	0.0	0.0	3	0.0	0.598487119987673	0.0	0.0	121.92490151987825;	%8 solar
	2	258941.26111912684	103576.50444765073	3	0.0	1321.1288832608511	0.0	0.0	339.1360592231248;	%9 CCGT
	2	337616.59274195385	135046.63709678155	3	0.0	1722.5336364385398	0.0	0.0	154.43644354494006;	%9 OCGT
	2	2284.2548411548487	913.7019364619393	3	0.0	447.1040988754841	0.0	0.0	952.339264917961;	%9 biomass
	2	0.0	0.0	3	0.0	0.7432292129472646	0.0	0.0	696.8395944681637;	%9 offwind
	2	0.0	0.0	3	0.0	0.7493500400124025	0.0	0.0	329.06126034235035;	%9 onwind
	2	0.0	0.0	3	0.0	0.5997416744131759	0.0	0.0	121.92490151987825;	%9 solar
	2	0.0	0.0	3	0.0	0.2900426077964185	0.0	0.0	5320356.498584592;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.21753195584731388	0.0	0.0	5320356.498584592;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3159958916828985	0.0	0.0	5320356.498584592;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.23699691876217388	0.0	0.0	5320356.498584592;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3297247681928945	0.0	0.0	5320356.498584592;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.24729357614467085	0.0	0.0	5320356.498584592;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.27535332286592107	0.0	0.0	5320356.498584592;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.2065149921494408	0.0	0.0	5320356.498584592;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3275906009733704	0.0	0.0	5320356.498584592;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.24569295073002778	0.0	0.0	5320356.498584592;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31514237987082666	0.0	0.0	5320356.498584592;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.23635678490312	0.0	0.0	5320356.498584592;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2999737134805074	0.0	0.0	5320356.498584592;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.22498028511038054	0.0	0.0	5320356.498584592;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3044957994645151	0.0	0.0	5320356.498584592;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.30106859754529913	0.0	0.0	5320356.498584592;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.22580144815897435	0.0	0.0	5320356.498584592;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2931730597008226	0.0	0.0	5320356.498584592;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.21987979477561692	0.0	0.0	5320356.498584592;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.28778444459603386	0.0	0.0	5320356.498584592;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.2158383334470254	0.0	0.0	5320356.498584592;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2755957587132579	0.0	0.0	5320356.498584592;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.20669681903494344	0.0	0.0	5320356.498584592;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.30397356252557967	0.0	0.0	5320356.498584592;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.22798017189418474	0.0	0.0	5320356.498584592;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.27903287771112056	0.0	0.0	5320356.498584592;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.2092746582833404	0.0	0.0	5320356.498584592;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.2807115967455776	0.0	0.0	5320356.498584592;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.21053369755918322	0.0	0.0	5320356.498584592;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3253090526527645	0.0	0.0	5320356.498584592;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.2439817894895734	0.0	0.0	5320356.498584592;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.297795189068741	0.0	0.0	5320356.498584592;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.22334639180155574	0.0	0.0	5320356.498584592;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3095007774591605	0.0	0.0	5320356.498584592;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.23212558309437037	0.0	0.0	5320356.498584592;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3137712682472393	0.0	0.0	5320356.498584592;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.23532845118542944	0.0	0.0	5320356.498584592;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.27859179927250693	0.0	0.0	5320356.498584592;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.32427789716122707	0.0	0.0	5320356.498584592;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.24320842287092032	0.0	0.0	5320356.498584592;	%DE1 76 PHS (pump mode)
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

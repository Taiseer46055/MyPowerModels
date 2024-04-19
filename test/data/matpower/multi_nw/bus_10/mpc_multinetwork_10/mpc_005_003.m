function mpc = mpc_005_003
%MPC_005_003	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 5 	Weight: 34
%	Time step: 3

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 34;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1464.14447129718	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	7775.846722960086	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4711.928447381428	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2991.06717972623	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	7412.1613915574235	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	9266.606848342955	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	2501.9991398872917	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4359.703743718252	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6036.332506361644	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	4410.955128660921	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.623418733404869	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.303872244786836	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.819552168542173	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.793239313065237	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.164412772047159	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.727523287185598	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.3757372232006793	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.03454457268254	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.671498369102697	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.112657279306637	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.929551270535689	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.640113372256097	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	16.662213673005432	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	283558.6324710038	113423.45298840152	3	0.0	1446.7277166887948	0.0	0.0	384.3542004528748;	%0 CCGT
	2	382525.47565981245	153010.19026392497	3	0.0	1951.6605901010835	0.0	0.0	175.02796935093204;	%0 OCGT
	2	2588.893075836218	1035.5572303344873	3	0.0	506.7318606060321	0.0	0.0	1079.317833573689;	%0 biomass
	2	0.0	0.0	3	0.0	0.8400813297836301	0.0	0.0	746.814757908944;	%0 offwind
	2	0.0	0.0	3	0.0	0.8356261195567983	0.0	0.0	372.9360950546637;	%0 onwind
	2	0.0	0.0	3	0.0	0.6686971025734494	0.0	0.0	138.18155505586202;	%0 solar
	2	284832.4131449443	113932.96525797772	3	0.0	1453.2265976782874	0.0	0.0	384.3542004528748;	%1 CCGT
	2	366351.7556164501	146540.70224658	3	0.0	1869.1416102880103	0.0	0.0	175.02796935093204;	%1 OCGT
	2	2588.838558888617	1035.5354235554469	3	0.0	506.72118983922826	0.0	0.0	1079.317833573689;	%1 biomass
	2	0.0	0.0	3	0.0	0.8504955792035117	0.0	0.0	372.9360950546637;	%1 onwind
	2	0.0	0.0	3	0.0	0.6844421116616622	0.0	0.0	138.18155505586202;	%1 solar
	2	271531.51572814974	108612.60629125989	3	0.0	1385.3648761640288	0.0	0.0	384.3542004528748;	%2 CCGT
	2	376422.16667309846	150568.8666692394	3	0.0	1920.5212585362165	0.0	0.0	175.02796935093204;	%2 OCGT
	2	2588.817211581205	1035.526884632482	3	0.0	506.7170114662761	0.0	0.0	1079.3178335736889;	%2 biomass
	2	0.0	0.0	3	0.0	0.8538425215141368	0.0	0.0	372.9360950546637;	%2 onwind
	2	0.0	0.0	3	0.0	0.6786671784654048	0.0	0.0	138.18155505586202;	%2 solar
	2	251426.22750790685	100570.49100316274	3	0.0	1282.7868750403409	0.0	0.0	384.3542004528748;	%3 CCGT
	2	361183.88530850504	144473.554123402	3	0.0	1842.7749250433926	0.0	0.0	175.02796935093204;	%3 OCGT
	2	2588.890915910465	1035.5563663641858	3	0.0	506.73143783724106	0.0	0.0	1079.317833573689;	%3 biomass
	2	0.0	0.0	3	0.0	0.8549255804487056	0.0	0.0	764.6741998147546;	%3 offwind
	2	0.0	0.0	3	0.0	0.8396266895724996	0.0	0.0	372.9360950546638;	%3 onwind
	2	0.0	0.0	3	0.0	0.6749107746112738	0.0	0.0	138.18155505586202;	%3 solar
	2	269223.86392463715	107689.54556985485	3	0.0	1373.5911424726382	0.0	0.0	384.3542004528748;	%4 CCGT
	2	371570.23660318437	148628.09464127375	3	0.0	1895.7665132815528	0.0	0.0	175.0279693509321;	%4 OCGT
	2	2588.893721642827	1035.5574886571308	3	0.0	506.73198701171003	0.0	0.0	1079.317833573689;	%4 biomass
	2	0.0	0.0	3	0.0	0.8418667999468696	0.0	0.0	372.9360950546637;	%4 onwind
	2	0.0	0.0	3	0.0	0.6648836353316124	0.0	0.0	138.18155505586202;	%4 solar
	2	280310.6667928772	112124.26671715088	3	0.0	1430.1564632289653	0.0	0.0	384.3542004528748;	%5 CCGT
	2	373293.0000955153	149317.20003820612	3	0.0	1904.5561229363025	0.0	0.0	175.02796935093204;	%5 OCGT
	2	2588.9061855602718	1035.5624742241089	3	0.0	506.73442661191467	0.0	0.0	1079.317833573689;	%5 biomass
	2	0.0	0.0	3	0.0	0.8505988914792149	0.0	0.0	372.9360950546637;	%5 onwind
	2	0.0	0.0	3	0.0	0.6882134053828558	0.0	0.0	138.18155505586202;	%5 solar
	2	266613.60900041024	106645.4436001641	3	0.0	1360.2735153082153	0.0	0.0	384.3542004528748;	%6 CCGT
	2	376599.4863205657	150639.79452822625	3	0.0	1921.4259506151307	0.0	0.0	175.0279693509321;	%6 OCGT
	2	2588.86993057649	1035.5479722305959	3	0.0	506.72733031444307	0.0	0.0	1079.317833573689;	%6 biomass
	2	0.0	0.0	3	0.0	0.8625269330099695	0.0	0.0	372.9360950546637;	%6 onwind
	2	0.0	0.0	3	0.0	0.6868579779032987	0.0	0.0	138.181555055862;	%6 solar
	2	270952.5383600299	108381.01534401196	3	0.0	1382.4109100001524	0.0	0.0	384.3542004528748;	%7 CCGT
	2	379112.2659295446	151644.9063718178	3	0.0	1934.2462547425741	0.0	0.0	175.02796935093204;	%7 OCGT
	2	2588.903861927173	1035.561544770869	3	0.0	506.7339718001904	0.0	0.0	1079.317833573689;	%7 biomass
	2	0.0	0.0	3	0.0	0.8558541853527175	0.0	0.0	372.9360950546637;	%7 onwind
	2	0.0	0.0	3	0.0	0.6747778850464267	0.0	0.0	138.18155505586202;	%7 solar
	2	291406.1613579287	116562.46454317149	3	0.0	1486.7661293771873	0.0	0.0	384.3542004528748;	%8 CCGT
	2	381197.4591505764	152478.98366023056	3	0.0	1944.8849956662057	0.0	0.0	175.02796935093204;	%8 OCGT
	2	2588.9240410502657	1035.5696164201063	3	0.0	506.73792152089754	0.0	0.0	1079.317833573689;	%8 biomass
	2	0.0	0.0	3	0.0	0.8616160407495415	0.0	0.0	372.9360950546637;	%8 onwind
	2	0.0	0.0	3	0.0	0.678285402652696	0.0	0.0	138.18155505586202;	%8 solar
	2	293466.7626016771	117386.70504067083	3	0.0	1497.2794010289645	0.0	0.0	384.3542004528748;	%9 CCGT
	2	382632.138440881	153052.8553763524	3	0.0	1952.2047879636784	0.0	0.0	175.0279693509321;	%9 OCGT
	2	2588.822153308828	1035.5288613235314	3	0.0	506.7179787255487	0.0	0.0	1079.317833573689;	%9 biomass
	2	0.0	0.0	3	0.0	0.8423264413402332	0.0	0.0	789.7515403972521;	%9 offwind
	2	0.0	0.0	3	0.0	0.8492633786807229	0.0	0.0	372.9360950546637;	%9 onwind
	2	0.0	0.0	3	0.0	0.6797072310015994	0.0	0.0	138.18155505586202;	%9 solar
	2	0.0	0.0	3	0.0	0.3287149555026076	0.0	0.0	6029737.365062538;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.24653621662695568	0.0	0.0	6029737.365062538;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3581286772406183	0.0	0.0	6029737.365062538;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.26859650793046375	0.0	0.0	6029737.365062538;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3736880706186137	0.0	0.0	6029737.365062538;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.2802660529639603	0.0	0.0	6029737.365062538;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3120670992480439	0.0	0.0	6029737.365062538;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.2340503244360329	0.0	0.0	6029737.365062538;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37126934776981974	0.0	0.0	6029737.365062538;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.2784520108273648	0.0	0.0	6029737.365062538;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3571613638536036	0.0	0.0	6029737.365062538;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.2678710228902027	0.0	0.0	6029737.365062538;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3399702086112417	0.0	0.0	6029737.365062538;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.2549776564584313	0.0	0.0	6029737.365062538;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.34509523939311715	0.0	0.0	6029737.365062538;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.3412110772180057	0.0	0.0	6029737.365062538;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.2559083079135043	0.0	0.0	6029737.365062538;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3322628009942656	0.0	0.0	6029737.365062538;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.24919710074569917	0.0	0.0	6029737.365062538;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.326155703875505	0.0	0.0	6029737.365062538;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.24461677790662878	0.0	0.0	6029737.365062538;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3123418598750256	0.0	0.0	6029737.365062538;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.2342563949062692	0.0	0.0	6029737.365062538;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3445033708623236	0.0	0.0	6029737.365062538;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.2583775281467427	0.0	0.0	6029737.365062538;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3162372614059366	0.0	0.0	6029737.365062538;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.23717794605445247	0.0	0.0	6029737.365062538;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.31813980964498795	0.0	0.0	6029737.365062538;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.23860485723374097	0.0	0.0	6029737.365062538;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3686835930064664	0.0	0.0	6029737.365062538;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.2765126947548498	0.0	0.0	6029737.365062538;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3375012142779064	0.0	0.0	6029737.365062538;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.2531259107084298	0.0	0.0	6029737.365062538;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.35076754778704855	0.0	0.0	6029737.365062538;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.26307566084028644	0.0	0.0	6029737.365062538;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3556074373468712	0.0	0.0	6029737.365062538;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.26670557801015343	0.0	0.0	6029737.365062538;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3157373725088412	0.0	0.0	6029737.365062538;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.36751495011605734	0.0	0.0	6029737.365062538;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.275636212587043	0.0	0.0	6029737.365062538;	%DE1 76 PHS (pump mode)
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

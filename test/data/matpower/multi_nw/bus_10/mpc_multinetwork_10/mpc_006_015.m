function mpc = mpc_006_015
%MPC_006_015	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 6 	Weight: 20
%	Time step: 15

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 20;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	747.5938145239273	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	9577.978586499721	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	6217.9570622168285	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2783.56527491027	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	6654.209909819276	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	7951.6720259865315	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	12520.063995045235	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	5750.802767965926	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5606.375327115718	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	11227.979829623488	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	12.76497197761797	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.085673884436483	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.14079169681704914	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6578724882364978	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.9219910537417644	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.9630194738229372	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6855359812829072	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.302283525386159	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	11.628381492060003	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.3440333054889311	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.6484400820996599	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.8172299888970908	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.1520730452972342	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	1.0958302605970554	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	2.8672339235543824	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.7351765225115334	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.931810723310572	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.31358328261491075	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.1820371049063	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.19285728987635975	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	13.320694149123181	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.60076443779762	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.24338564446895838	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	166799.19557117872	66719.67822847149	3	0.0	851.0163039345853	0.0	0.0	50644.318177319976;	%0 CCGT
	2	225014.9856822426	90005.99427289705	3	0.0	1148.0356412359315	0.0	0.0	23062.50890271105;	%0 OCGT
	2	1522.8782799036576	609.1513119614631	3	0.0	298.07756506237183	0.0	0.0	6348.928432786406;	%0 biomass
	2	0.0	0.0	3	0.0	0.4941654881080177	0.0	0.0	7468.147579089441;	%0 offwind
	2	0.0	0.0	3	0.0	0.49154477620988135	0.0	0.0	8336.21859533954;	%0 onwind
	2	0.0	0.0	3	0.0	0.3933512368079114	0.0	0.0	2275.9314950377275;	%0 solar
	2	167548.47832055547	67019.39132822219	3	0.0	854.8391751048748	0.0	0.0	50644.318177319976;	%1 CCGT
	2	215501.03271555886	86200.41308622355	3	0.0	1099.4950648753002	0.0	0.0	23062.50890271105;	%1 OCGT
	2	1522.8462111109513	609.1384844443805	3	0.0	298.0712881407225	0.0	0.0	6348.928432786406;	%1 biomass
	2	0.0	0.0	3	0.0	0.5002915171785363	0.0	0.0	8336.21859533954;	%1 onwind
	2	0.0	0.0	3	0.0	0.40261300685980134	0.0	0.0	2275.9314950377275;	%1 solar
	2	159724.42101655866	63889.768406623465	3	0.0	814.9205153906053	0.0	0.0	50644.318177319976;	%2 CCGT
	2	221424.80392535203	88569.92157014081	3	0.0	1129.718387374245	0.0	0.0	23062.50890271105;	%2 OCGT
	2	1522.833653871297	609.1334615485188	3	0.0	298.06883027428006	0.0	0.0	6348.928432786406;	%2 biomass
	2	0.0	0.0	3	0.0	0.5022603067730216	0.0	0.0	8336.21859533954;	%2 onwind
	2	0.0	0.0	3	0.0	0.3992159873325911	0.0	0.0	2275.9314950377275;	%2 solar
	2	147897.78088700402	59159.11235480161	3	0.0	754.5805147296123	0.0	0.0	50644.318177319976;	%3 CCGT
	2	212461.10900500295	84984.44360200118	3	0.0	1083.9852500255251	0.0	0.0	23062.50890271105;	%3 OCGT
	2	1522.877009359097	609.1508037436388	3	0.0	298.07731637484767	0.0	0.0	6348.928432786406;	%3 biomass
	2	0.0	0.0	3	0.0	0.5028974002639445	0.0	0.0	7646.741998147546;	%3 offwind
	2	0.0	0.0	3	0.0	0.4938980526897057	0.0	0.0	8336.218595339542;	%3 onwind
	2	0.0	0.0	3	0.0	0.39700633800663165	0.0	0.0	2275.9314950377275;	%3 solar
	2	158366.9787791983	63346.79151167932	3	0.0	807.9947896897872	0.0	0.0	50644.318177319976;	%4 CCGT
	2	218570.72741363785	87428.29096545515	3	0.0	1115.1567725185605	0.0	0.0	23062.50890271105;	%4 OCGT
	2	1522.8786597898984	609.1514639159593	3	0.0	298.07763941865295	0.0	0.0	6348.928432786406;	%4 biomass
	2	0.0	0.0	3	0.0	0.4952157646746292	0.0	0.0	8336.21859533954;	%4 onwind
	2	0.0	0.0	3	0.0	0.3911080207833014	0.0	0.0	2275.9314950377275;	%4 solar
	2	164888.6275252219	65955.45101008876	3	0.0	841.2685077817443	0.0	0.0	50644.318177319976;	%5 CCGT
	2	219584.11770324426	87833.6470812977	3	0.0	1120.3271311390015	0.0	0.0	23062.50890271105;	%5 OCGT
	2	1522.8859915060423	609.1543966024169	3	0.0	298.07907447759686	0.0	0.0	6348.928432786406;	%5 biomass
	2	0.0	0.0	3	0.0	0.5003522891054205	0.0	0.0	8336.21859533954;	%5 onwind
	2	0.0	0.0	3	0.0	0.40483141493109165	0.0	0.0	2275.9314950377275;	%5 solar
	2	156831.53470612367	62732.61388244947	3	0.0	800.1608913577738	0.0	0.0	50644.318177319976;	%6 CCGT
	2	221529.10960033274	88611.6438401331	3	0.0	1130.250559185371	0.0	0.0	23062.50890271105;	%6 OCGT
	2	1522.8646650449941	609.1458660179976	3	0.0	298.0749001849665	0.0	0.0	6348.928432786406;	%6 biomass
	2	0.0	0.0	3	0.0	0.5073687841235115	0.0	0.0	8336.21859533954;	%6 onwind
	2	0.0	0.0	3	0.0	0.40403410464899925	0.0	0.0	2275.931495037727;	%6 solar
	2	159383.84609413525	63753.5384376541	3	0.0	813.1828882353839	0.0	0.0	50644.318177319976;	%7 CCGT
	2	223007.21525267325	89202.8861010693	3	0.0	1137.7919145544554	0.0	0.0	23062.50890271105;	%7 OCGT
	2	1522.8846246630428	609.1538498652171	3	0.0	298.07880694128846	0.0	0.0	6348.928432786406;	%7 biomass
	2	0.0	0.0	3	0.0	0.503443638442775	0.0	0.0	8336.21859533954;	%7 onwind
	2	0.0	0.0	3	0.0	0.3969281676743686	0.0	0.0	2275.9314950377275;	%7 solar
	2	171415.38903407572	68566.15561363028	3	0.0	874.5683113983455	0.0	0.0	50644.318177319976;	%8 CCGT
	2	224233.79950033905	89693.51980013562	3	0.0	1144.0499974507093	0.0	0.0	23062.50890271105;	%8 OCGT
	2	1522.8964947354505	609.1585978941802	3	0.0	298.08113030641033	0.0	0.0	6348.928432786406;	%8 biomass
	2	0.0	0.0	3	0.0	0.5068329651467891	0.0	0.0	8336.21859533954;	%8 onwind
	2	0.0	0.0	3	0.0	0.39899141332511534	0.0	0.0	2275.9314950377275;	%8 solar
	2	172627.5074127512	69051.00296510049	3	0.0	880.7525888405673	0.0	0.0	50644.318177319976;	%9 CCGT
	2	225077.7284946359	90031.09139785435	3	0.0	1148.3557576256933	0.0	0.0	23062.50890271105;	%9 OCGT
	2	1522.836560769899	609.1346243079596	3	0.0	298.06939925032276	0.0	0.0	6348.928432786406;	%9 biomass
	2	0.0	0.0	3	0.0	0.4954861419648431	0.0	0.0	7897.515403972521;	%9 offwind
	2	0.0	0.0	3	0.0	0.49956669334160164	0.0	0.0	8336.21859533954;	%9 onwind
	2	0.0	0.0	3	0.0	0.3998277829421173	0.0	0.0	2275.9314950377275;	%9 solar
	2	0.0	0.0	3	0.0	0.19336173853094565	0.0	0.0	404.8977548390101;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.14502130389820925	0.0	0.0	404.8977548390101;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.210663927788599	0.0	0.0	404.8977548390101;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.15799794584144924	0.0	0.0	404.8977548390101;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21981651212859632	0.0	0.0	404.8977548390101;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.16486238409644724	0.0	0.0	404.8977548390101;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18356888191061405	0.0	0.0	404.8977548390101;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.13767666143296053	0.0	0.0	404.8977548390101;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21839373398224693	0.0	0.0	404.8977548390101;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.1637953004866852	0.0	0.0	404.8977548390101;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21009491991388446	0.0	0.0	404.8977548390101;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.15757118993541336	0.0	0.0	404.8977548390101;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1999824756536716	0.0	0.0	404.8977548390101;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.1499868567402537	0.0	0.0	404.8977548390101;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20299719964301008	0.0	0.0	404.8977548390101;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.20071239836353277	0.0	0.0	404.8977548390101;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.15053429877264957	0.0	0.0	404.8977548390101;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19544870646721504	0.0	0.0	404.8977548390101;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.1465865298504113	0.0	0.0	404.8977548390101;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1918562963973559	0.0	0.0	404.8977548390101;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.14389222229801693	0.0	0.0	404.8977548390101;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1837305058088386	0.0	0.0	404.8977548390101;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.13779787935662896	0.0	0.0	404.8977548390101;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20264904168371978	0.0	0.0	404.8977548390101;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.15198678126278983	0.0	0.0	404.8977548390101;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18602191847408037	0.0	0.0	404.8977548390101;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.13951643885556028	0.0	0.0	404.8977548390101;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.18714106449705173	0.0	0.0	404.8977548390101;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.1403557983727888	0.0	0.0	404.8977548390101;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.21687270176850965	0.0	0.0	404.8977548390101;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.16265452632638222	0.0	0.0	404.8977548390101;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.19853012604582732	0.0	0.0	404.8977548390101;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.1488975945343705	0.0	0.0	404.8977548390101;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20633385163944032	0.0	0.0	404.8977548390101;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.15475038872958025	0.0	0.0	404.8977548390101;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.20918084549815952	0.0	0.0	404.8977548390101;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.15688563412361964	0.0	0.0	404.8977548390101;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.1857278661816713	0.0	0.0	404.8977548390101;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.21618526477415137	0.0	0.0	404.8977548390101;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.16213894858061353	0.0	0.0	404.8977548390101;	%DE1 76 PHS (pump mode)
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

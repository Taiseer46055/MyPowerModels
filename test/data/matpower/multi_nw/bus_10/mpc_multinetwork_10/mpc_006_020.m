function mpc = mpc_006_020
%MPC_006_020	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 6 	Weight: 47
%	Time step: 20

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 47;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	940.4210767323897	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	6111.247028574406	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	5729.0422858206	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	2629.52855939094	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	2668.103952035713	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	443.0174728575865	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	5526.980933055577	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4123.310955926994	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	6383.117970166457	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	5173.7950973627385	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.8900271144646	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.132647237089469	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142859	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.991835787736166	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.930004177215192	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.460896291279643	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	7.667898942869358	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.395387182954048	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.512424130066765	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.949334269662922	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.494125134172864	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49474999999998	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.019351984494808	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.298963432313621	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.047139226743301	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538467	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.492832005960958	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.188815522355116	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	391978.10959227	156791.243836908	3	0.0	1999.8883142462753	0.0	0.0	531.3131594495622;	%0 CCGT
	2	528785.2163532701	211514.08654130806	3	0.0	2697.883756904439	0.0	0.0	241.95042822040608;	%0 OCGT
	2	3578.7639577735954	1431.5055831094382	3	0.0	700.4822778965738	0.0	0.0	1491.9981817048056;	%0 biomass
	2	0.0	0.0	3	0.0	1.1612888970538415	0.0	0.0	1032.3615771094226;	%0 offwind
	2	0.0	0.0	3	0.0	1.1551302240932213	0.0	0.0	515.5293078696822;	%0 onwind
	2	0.0	0.0	3	0.0	0.9243754064985917	0.0	0.0	191.01567904780924;	%0 solar
	2	393738.92405330535	157495.56962132215	3	0.0	2008.872061496456	0.0	0.0	531.3131594495622;	%1 CCGT
	2	506427.4268815633	202570.97075262532	3	0.0	2583.8134024569554	0.0	0.0	241.95042822040608;	%1 OCGT
	2	3578.6885961107355	1431.4754384442942	3	0.0	700.4675271306979	0.0	0.0	1491.9981817048056;	%1 biomass
	2	0.0	0.0	3	0.0	1.1756850653695603	0.0	0.0	515.5293078696822;	%1 onwind
	2	0.0	0.0	3	0.0	0.9461405661205331	0.0	0.0	191.01567904780924;	%1 solar
	2	375352.3893889129	150140.95575556514	3	0.0	1915.0632111679224	0.0	0.0	531.3131594495622;	%2 CCGT
	2	520348.28922457725	208139.31568983092	3	0.0	2654.8382103294757	0.0	0.0	241.95042822040608;	%2 OCGT
	2	3578.6590865975477	1431.4636346390191	3	0.0	700.4617511445582	0.0	0.0	1491.9981817048053;	%2 biomass
	2	0.0	0.0	3	0.0	1.1803117209166007	0.0	0.0	515.5293078696822;	%2 onwind
	2	0.0	0.0	3	0.0	0.9381575702315891	0.0	0.0	191.01567904780924;	%2 solar
	2	347559.78508445946	139023.91403378377	3	0.0	1773.264209614589	0.0	0.0	531.3131594495622;	%3 CCGT
	2	499283.6061617569	199713.44246470276	3	0.0	2547.365337559984	0.0	0.0	241.95042822040608;	%3 OCGT
	2	3578.760971993878	1431.504388797551	3	0.0	700.4816934808921	0.0	0.0	1491.9981817048056;	%3 biomass
	2	0.0	0.0	3	0.0	1.1818088906202695	0.0	0.0	1057.0496291556904;	%3 offwind
	2	0.0	0.0	3	0.0	1.1606604238208085	0.0	0.0	515.5293078696823;	%3 onwind
	2	0.0	0.0	3	0.0	0.9329648943155844	0.0	0.0	191.01567904780924;	%3 solar
	2	372162.40013111604	148864.9600524464	3	0.0	1898.7877557709999	0.0	0.0	531.3131594495622;	%4 CCGT
	2	513641.209422049	205456.4837688196	3	0.0	2620.618415418617	0.0	0.0	241.9504282204061;	%4 OCGT
	2	3578.764850506261	1431.5059402025042	3	0.0	700.4824526338344	0.0	0.0	1491.9981817048056;	%4 biomass
	2	0.0	0.0	3	0.0	1.1637570469853786	0.0	0.0	515.5293078696822;	%4 onwind
	2	0.0	0.0	3	0.0	0.9191038488407584	0.0	0.0	191.01567904780924;	%4 solar
	2	387488.2746842714	154995.30987370858	3	0.0	1976.9809932870992	0.0	0.0	531.3131594495622;	%5 CCGT
	2	516022.676602624	206409.07064104962	3	0.0	2632.7687581766536	0.0	0.0	241.95042822040608;	%5 OCGT
	2	3578.782080039199	1431.5128320156798	3	0.0	700.4858250223526	0.0	0.0	1491.9981817048056;	%5 biomass
	2	0.0	0.0	3	0.0	1.1758278793977381	0.0	0.0	515.5293078696822;	%5 onwind
	2	0.0	0.0	3	0.0	0.9513538250880653	0.0	0.0	191.01567904780924;	%5 solar
	2	368554.1065593906	147421.64262375626	3	0.0	1880.3780946907684	0.0	0.0	531.3131594495622;	%6 CCGT
	2	520593.4075607819	208237.36302431277	3	0.0	2656.088814085622	0.0	0.0	241.9504282204061;	%6 OCGT
	2	3578.731962855736	1431.4927851422942	3	0.0	700.4760154346714	0.0	0.0	1491.9981817048056;	%6 biomass
	2	0.0	0.0	3	0.0	1.1923166426902518	0.0	0.0	515.5293078696822;	%6 onwind
	2	0.0	0.0	3	0.0	0.9494801459251483	0.0	0.0	191.01567904780921;	%6 solar
	2	374552.0383212178	149820.81532848714	3	0.0	1910.979787353152	0.0	0.0	531.3131594495622;	%7 CCGT
	2	524066.95584378217	209626.78233751288	3	0.0	2673.8109992029704	0.0	0.0	241.95042822040608;	%7 OCGT
	2	3578.778867958151	1431.51154718326	3	0.0	700.4851963120278	0.0	0.0	1491.9981817048056;	%7 biomass
	2	0.0	0.0	3	0.0	1.1830925503405212	0.0	0.0	515.5293078696822;	%7 onwind
	2	0.0	0.0	3	0.0	0.9327811940347662	0.0	0.0	191.01567904780924;	%7 solar
	2	402826.16423007793	161130.46569203117	3	0.0	2055.235531786112	0.0	0.0	531.3131594495622;	%8 CCGT
	2	526949.4288257968	210779.7715303187	3	0.0	2688.5174940091665	0.0	0.0	241.95042822040608;	%8 OCGT
	2	3578.8067626283087	1431.5227050513233	3	0.0	700.4906562200642	0.0	0.0	1491.9981817048056;	%8 biomass
	2	0.0	0.0	3	0.0	1.1910574680949544	0.0	0.0	515.5293078696822;	%8 onwind
	2	0.0	0.0	3	0.0	0.937629821314021	0.0	0.0	191.01567904780924;	%8 solar
	2	405674.64241996536	162269.85696798615	3	0.0	2069.7685837753334	0.0	0.0	531.3131594495622;	%9 CCGT
	2	528932.6619623944	211573.06478495774	3	0.0	2698.636030420379	0.0	0.0	241.9504282204061;	%9 OCGT
	2	3578.6659178092627	1431.4663671237051	3	0.0	700.4630882382585	0.0	0.0	1491.9981817048056;	%9 biomass
	2	0.0	0.0	3	0.0	1.1643924336173812	0.0	0.0	1091.7153646667898;	%9 offwind
	2	0.0	0.0	3	0.0	1.173981729352764	0.0	0.0	515.5293078696822;	%9 onwind
	2	0.0	0.0	3	0.0	0.9395952899139756	0.0	0.0	191.01567904780924;	%9 solar
	2	0.0	0.0	3	0.0	0.45440008554772227	0.0	0.0	8335225.181115861;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.34080006416079167	0.0	0.0	8335225.181115861;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.49506023030320767	0.0	0.0	8335225.181115861;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.37129517272740575	0.0	0.0	8335225.181115861;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5165688035022014	0.0	0.0	8335225.181115861;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.387426602626651	0.0	0.0	8335225.181115861;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.431386872489943	0.0	0.0	8335225.181115861;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.32354015436745726	0.0	0.0	8335225.181115861;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5132252748582803	0.0	0.0	8335225.181115861;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.3849189561437102	0.0	0.0	8335225.181115861;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4937230617976285	0.0	0.0	8335225.181115861;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.37029229634822136	0.0	0.0	8335225.181115861;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46995881778612825	0.0	0.0	8335225.181115861;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.3524691133395962	0.0	0.0	8335225.181115861;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4770434191610737	0.0	0.0	8335225.181115861;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.471674136154302	0.0	0.0	8335225.181115861;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.3537556021157265	0.0	0.0	8335225.181115861;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.45930446019795534	0.0	0.0	8335225.181115861;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.34447834514846654	0.0	0.0	8335225.181115861;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.45086229653378634	0.0	0.0	8335225.181115861;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.33814672240033977	0.0	0.0	8335225.181115861;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4317666886507707	0.0	0.0	8335225.181115861;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.32382501648807804	0.0	0.0	8335225.181115861;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4762252479567414	0.0	0.0	8335225.181115861;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.35716893596755606	0.0	0.0	8335225.181115861;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43715150841408884	0.0	0.0	8335225.181115861;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.3278636313105666	0.0	0.0	8335225.181115861;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43978150156807155	0.0	0.0	8335225.181115861;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.3298361261760537	0.0	0.0	8335225.181115861;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.5096508491559977	0.0	0.0	8335225.181115861;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.38223813686699826	0.0	0.0	8335225.181115861;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46654579620769415	0.0	0.0	8335225.181115861;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.3499093471557706	0.0	0.0	8335225.181115861;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.48488455135268477	0.0	0.0	8335225.181115861;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.3636634135145136	0.0	0.0	8335225.181115861;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4915749869206749	0.0	0.0	8335225.181115861;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.3686812401905062	0.0	0.0	8335225.181115861;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.43646048552692757	0.0	0.0	8335225.181115861;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.5080353722192558	0.0	0.0	8335225.181115861;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.3810265291644418	0.0	0.0	8335225.181115861;	%DE1 76 PHS (pump mode)
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

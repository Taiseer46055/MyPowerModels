function mpc = mpc_010_009
%MPC_010_009	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 10 	Weight: 36
%	Time step: 9

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 36;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	2662.2892062755104	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	12635.284088923374	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	9500.899839641655	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	7013.825677454592	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	5268.269543077674	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	5949.524568541816	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	15034.248892710451	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	7722.700219894207	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	4340.925490857316	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	11375.866514453766	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0859453848135806	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.914011698690512	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.068229078059304	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666667	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.93000417721519	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.347583923124404	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.11	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.12585866747405694	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.423509654213486	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	19.03134203502625	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.8	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.627262599278453	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	17.448688009541062	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.425065824509822	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.029751484476830056	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.223241110968374	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.005803818276707249	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	18.189174489846003	0.0	0.0	221	221	2191	28.0	5;	%9 solar
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
	2	300238.55202812166	120095.42081124867	3	0.0	1531.8293470822534	0.0	0.0	91159.77271917595;	%0 CCGT
	2	405026.9742280367	162010.78969121468	3	0.0	2066.464154224677	0.0	0.0	41512.51602487989;	%0 OCGT
	2	2741.1809038265837	1096.4723615306336	3	0.0	536.5396171122693	0.0	0.0	11428.071179015531;	%0 biomass
	2	0.0	0.0	3	0.0	0.8894978785944319	0.0	0.0	13442.665642360993;	%0 offwind
	2	0.0	0.0	3	0.0	0.8847805971777865	0.0	0.0	15005.193471611174;	%0 onwind
	2	0.0	0.0	3	0.0	0.7080322262542404	0.0	0.0	4096.676691067909;	%0 solar
	2	301587.2609769999	120634.90439079994	3	0.0	1538.7105151887747	0.0	0.0	91159.77271917595;	%1 CCGT
	2	387901.85888800595	155160.74355520238	3	0.0	1979.0911167755403	0.0	0.0	41512.51602487989;	%1 OCGT
	2	2741.123179999712	1096.449271999885	3	0.0	536.5283186533005	0.0	0.0	11428.071179015531;	%1 biomass
	2	0.0	0.0	3	0.0	0.9005247309213653	0.0	0.0	15005.193471611174;	%1 onwind
	2	0.0	0.0	3	0.0	0.7247034123476424	0.0	0.0	4096.676691067909;	%1 solar
	2	287503.9578298056	115001.58313192223	3	0.0	1466.8569277030895	0.0	0.0	91159.77271917595;	%2 CCGT
	2	398564.6470656337	159425.85882625348	3	0.0	2033.4930972736408	0.0	0.0	41512.51602487989;	%2 OCGT
	2	2741.1005769683343	1096.4402307873338	3	0.0	536.5238944937041	0.0	0.0	11428.071179015531;	%2 biomass
	2	0.0	0.0	3	0.0	0.9040685521914389	0.0	0.0	15005.193471611174;	%2 onwind
	2	0.0	0.0	3	0.0	0.718588777198664	0.0	0.0	4096.676691067909;	%2 solar
	2	266216.0055966072	106486.4022386429	3	0.0	1358.244926513302	0.0	0.0	91159.77271917595;	%3 CCGT
	2	382429.9962090053	152971.9984836021	3	0.0	1951.1734500459452	0.0	0.0	41512.51602487989;	%3 OCGT
	2	2741.1786168463746	1096.4714467385497	3	0.0	536.5391694747258	0.0	0.0	11428.071179015531;	%3 biomass
	2	0.0	0.0	3	0.0	0.9052153204751001	0.0	0.0	13764.135596665583;	%3 offwind
	2	0.0	0.0	3	0.0	0.8890164948414703	0.0	0.0	15005.193471611175;	%3 onwind
	2	0.0	0.0	3	0.0	0.714611408411937	0.0	0.0	4096.676691067909;	%3 solar
	2	285060.561802557	114024.22472102278	3	0.0	1454.390621441617	0.0	0.0	91159.77271917595;	%4 CCGT
	2	393427.30934454815	157370.92373781925	3	0.0	2007.2821905334088	0.0	0.0	41512.51602487989;	%4 OCGT
	2	2741.181587621817	1096.4726350487267	3	0.0	536.5397509535753	0.0	0.0	11428.071179015531;	%4 biomass
	2	0.0	0.0	3	0.0	0.8913883764143326	0.0	0.0	15005.193471611174;	%4 onwind
	2	0.0	0.0	3	0.0	0.7039944374099425	0.0	0.0	4096.676691067909;	%4 solar
	2	296799.5295453994	118719.81181815975	3	0.0	1514.2833140071398	0.0	0.0	91159.77271917595;	%5 CCGT
	2	395251.4118658397	158100.56474633588	3	0.0	2016.5888360502026	0.0	0.0	41512.51602487989;	%5 OCGT
	2	2741.194784710876	1096.4779138843505	3	0.0	536.5423340596743	0.0	0.0	11428.071179015531;	%5 biomass
	2	0.0	0.0	3	0.0	0.900634120389757	0.0	0.0	15005.193471611174;	%5 onwind
	2	0.0	0.0	3	0.0	0.728696546875965	0.0	0.0	4096.676691067909;	%5 solar
	2	282296.7624710226	112918.70498840904	3	0.0	1440.2896044439929	0.0	0.0	91159.77271917595;	%6 CCGT
	2	398752.39728059893	159500.95891223956	3	0.0	2034.4510065336676	0.0	0.0	41512.51602487989;	%6 OCGT
	2	2741.1563970809893	1096.4625588323956	3	0.0	536.5348203329397	0.0	0.0	11428.071179015531;	%6 biomass
	2	0.0	0.0	3	0.0	0.9132638114223206	0.0	0.0	15005.193471611174;	%6 onwind
	2	0.0	0.0	3	0.0	0.7272613883681986	0.0	0.0	4096.676691067909;	%6 solar
	2	286890.92296944343	114756.36918777737	3	0.0	1463.729198823691	0.0	0.0	91159.77271917595;	%7 CCGT
	2	401412.98745481187	160565.19498192475	3	0.0	2048.02544619802	0.0	0.0	41512.51602487989;	%7 OCGT
	2	2741.192324393477	1096.4769297573907	3	0.0	536.5418524943192	0.0	0.0	11428.071179015531;	%7 biomass
	2	0.0	0.0	3	0.0	0.9061985491969949	0.0	0.0	15005.193471611174;	%7 onwind
	2	0.0	0.0	3	0.0	0.7144707018138635	0.0	0.0	4096.676691067909;	%7 solar
	2	308547.7002613363	123419.08010453451	3	0.0	1574.2229605170219	0.0	0.0	91159.77271917595;	%8 CCGT
	2	403620.8391006103	161448.33564024413	3	0.0	2059.289995411277	0.0	0.0	41512.51602487989;	%8 OCGT
	2	2741.213690523811	1096.4854762095242	3	0.0	536.5460345515386	0.0	0.0	11428.071179015531;	%8 biomass
	2	0.0	0.0	3	0.0	0.9122993372642203	0.0	0.0	15005.193471611174;	%8 onwind
	2	0.0	0.0	3	0.0	0.7181845439852076	0.0	0.0	4096.676691067909;	%8 solar
	2	310729.5133429522	124291.80533718088	3	0.0	1585.3546599130214	0.0	0.0	91159.77271917595;	%9 CCGT
	2	405139.9112903446	162055.96451613784	3	0.0	2067.040363726248	0.0	0.0	41512.51602487989;	%9 OCGT
	2	2741.1058093858182	1096.4423237543274	3	0.0	536.524918650581	0.0	0.0	11428.071179015531;	%9 biomass
	2	0.0	0.0	3	0.0	0.8918750555367175	0.0	0.0	14215.527727150538;	%9 offwind
	2	0.0	0.0	3	0.0	0.899220048014883	0.0	0.0	15005.193471611174;	%9 onwind
	2	0.0	0.0	3	0.0	0.7196900092958111	0.0	0.0	4096.676691067909;	%9 solar
	2	0.0	0.0	3	0.0	0.34805112935570215	0.0	0.0	728.8159587102181;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	-0.26103834701677664	0.0	0.0	728.8159587102181;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3791950700194782	0.0	0.0	728.8159587102181;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	-0.2843963025146087	0.0	0.0	728.8159587102181;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39566972183147336	0.0	0.0	728.8159587102181;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	-0.29675229137360504	0.0	0.0	728.8159587102181;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3304239874391053	0.0	0.0	728.8159587102181;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	-0.24781799057932896	0.0	0.0	728.8159587102181;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.39310872116804446	0.0	0.0	728.8159587102181;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	-0.29483154087603336	0.0	0.0	728.8159587102181;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.378170855844992	0.0	0.0	728.8159587102181;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	-0.28362814188374397	0.0	0.0	728.8159587102181;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3599684561766089	0.0	0.0	728.8159587102181;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	-0.2699763421324567	0.0	0.0	728.8159587102181;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.36539495935741817	0.0	0.0	728.8159587102181;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.361282317054359	0.0	0.0	728.8159587102181;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	-0.27096173779076926	0.0	0.0	728.8159587102181;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3518076716409871	0.0	0.0	728.8159587102181;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	-0.26385575373074033	0.0	0.0	728.8159587102181;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3453413335152406	0.0	0.0	728.8159587102181;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	-0.25900600013643044	0.0	0.0	728.8159587102181;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3307149104559095	0.0	0.0	728.8159587102181;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	-0.24803618284193213	0.0	0.0	728.8159587102181;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3647682750306956	0.0	0.0	728.8159587102181;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	-0.2735762062730217	0.0	0.0	728.8159587102181;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33483945325334463	0.0	0.0	728.8159587102181;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	-0.25112958994000845	0.0	0.0	728.8159587102181;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3368539160946931	0.0	0.0	728.8159587102181;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	-0.25264043707101985	0.0	0.0	728.8159587102181;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3903708631833174	0.0	0.0	728.8159587102181;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	-0.29277814738748803	0.0	0.0	728.8159587102181;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.35735422688248913	0.0	0.0	728.8159587102181;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	-0.26801567016186684	0.0	0.0	728.8159587102181;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3714009329509926	0.0	0.0	728.8159587102181;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	-0.2785506997132444	0.0	0.0	728.8159587102181;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.37652552189668714	0.0	0.0	728.8159587102181;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	-0.28239414142251534	0.0	0.0	728.8159587102181;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.33431015912700834	0.0	0.0	728.8159587102181;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.38913347659347247	0.0	0.0	728.8159587102181;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	-0.29185010744510437	0.0	0.0	728.8159587102181;	%DE1 76 PHS (pump mode)
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

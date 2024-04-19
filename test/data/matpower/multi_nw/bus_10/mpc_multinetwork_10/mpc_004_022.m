function mpc = mpc_004_022
%MPC_004_022	Extended MATPOWER case file generated from aggregated PyPSA generator and load data
%	Period: 4 	Weight: 43
%	Time step: 22

%%	MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 100;

%% Weight
mpc.weight = 43;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	1110.9456900023379	0.0	0.0	0.0	1	1.0	0.0	230.0	1	1.1	0.9;
	2	2	7960.536747207505	0.0	0.0	0.0	2	1.0	0.0	230.0	1	1.1	0.9;
	3	2	4886.649400779357	0.0	0.0	0.0	3	1.0	0.0	230.0	1	1.1	0.9;
	4	2	3024.3331706355734	0.0	0.0	0.0	4	1.0	0.0	230.0	1	1.1	0.9;
	5	2	6210.936217119917	0.0	0.0	0.0	5	1.0	0.0	230.0	1	1.1	0.9;
	6	2	7523.251822911347	0.0	0.0	0.0	6	1.0	0.0	230.0	1	1.1	0.9;
	7	2	6464.5567267676515	0.0	0.0	0.0	7	1.0	0.0	230.0	1	1.1	0.9;
	8	2	4833.213314560585	0.0	0.0	0.0	8	1.0	0.0	230.0	1	1.1	0.9;
	9	2	5343.9591287719895	0.0	0.0	0.0	9	1.0	0.0	230.0	1	1.1	0.9;
	10	2	6012.0486395679445	0.0	0.0	0.0	10	1.0	0.0	230.0	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	H	n0	nE	nE_max	P_b_nom	carrier
mpc.gen = [
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	183.0	78.4	5.0	1	1	91	224.0	0;	%0 CCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	135.7	78.4	2.8	1	1	91	224.0	1;	%0 OCGT
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.812136046511629	2.0436	2.8	43	43	2035	10.0	2;	%0 biomass
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.67520755664969	0.0	0.0	55	55	1198	17.0	3;	%0 offwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.76859903570107	0.0	0.0	146	146	536	38.0	4;	%0 onwind
	1	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	107	107	727	28.0	5;	%0 solar
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	218.276	78.4	5.0	10	10	110	224.0	0;	%1 CCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	151.275	78.4	2.8	3	3	110	224.0	1;	%1 OCGT
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.809952857142857	2.0436	2.8	42	42	2453	10.0	2;	%1 biomass
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.90379149373623	0.0	0.0	126	126	646	38.0	4;	%1 onwind
	2	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	160	160	876	28.0	5;	%1 solar
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	157.36666666666665	78.4	5.0	3	3	93	224.0	0;	%2 CCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	136.775	78.4	2.8	2	2	93	224.0	1;	%2 OCGT
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.930004177215192	2.0436	2.8	79	79	2067	10.0	2;	%2 biomass
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.23256085017743	0.0	0.0	79	79	544	38.0	4;	%2 onwind
	3	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	209	209	739	28.0	5;	%2 solar
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.0	78.4	5.0	2	2	138	224.0	0;	%3 CCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	120.10999999999999	78.4	2.8	2	2	138	224.0	1;	%3 OCGT
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.931164764150944	2.0436	2.8	106	106	3084	10.0	2;	%3 biomass
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	6.463147696972904	0.0	0.0	53	53	1814	17.0	3;	%3 offwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.734693679018074	0.0	0.0	257	257	812	38.0	4;	%3 onwind
	4	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	112	112	1102	28.0	5;	%3 solar
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	208.76153846153846	78.4	5.0	13	13	127	224.0	0;	%4 CCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	187.555	78.4	2.8	2	2	127	224.0	1;	%4 OCGT
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.925124485981307	2.0436	2.8	107	107	2827	10.0	2;	%4 biomass
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	3.9812772014122495	0.0	0.0	33	33	744	38.0	4;	%4 onwind
	5	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	313	313	1010	28.0	5;	%4 solar
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	195.79999999999998	78.4	5.0	2	2	103	224.0	0;	%5 CCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	181.14283333333333	78.4	2.8	3	3	103	224.0	1;	%5 OCGT
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.94933426966292	2.0436	2.8	89	89	2298	10.0	2;	%5 biomass
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.149647038128663	0.0	0.0	61	61	605	38.0	4;	%5 onwind
	6	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	264	264	821	28.0	5;	%5 solar
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	217.46415789473681	78.4	5.0	19	19	105	224.0	0;	%6 CCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.49475	78.4	2.8	4	4	105	224.0	1;	%6 OCGT
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.893259292682927	2.0436	2.8	41	41	2351	10.0	2;	%6 biomass
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	10.499409805941555	0.0	0.0	84	84	619	38.0	4;	%6 onwind
	7	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	118	118	840	28.0	5;	%6 solar
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	191.66000000000003	78.4	5.0	5	5	126	224.0	0;	%7 CCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	163.96	78.4	2.8	2	2	126	224.0	1;	%7 OCGT
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.906306436781607	2.0436	2.8	87	87	2805	10.0	2;	%7 biomass
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	4.9180241944526895	0.0	0.0	188	188	738	38.0	4;	%7 onwind
	8	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	170	170	1002	28.0	5;	%7 solar
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	193.0157142857143	78.4	5.0	7	7	144	224.0	0;	%8 CCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	169.005	78.4	2.8	4	4	144	224.0	1;	%8 OCGT
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.908332467532468	2.0436	2.8	77	77	3215	10.0	2;	%8 biomass
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	5.583490252989092	0.0	0.0	166	166	846	38.0	4;	%8 onwind
	9	0.0	0.0	0.0	-0.0	1.0	100.0	1	0.0	0.0	0.0	250	250	1149	28.0	5;	%8 solar
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	220.05538461538464	78.4	5.0	26	26	274	224.0	0;	%9 CCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	203.98000000000002	78.4	2.8	1	1	274	224.0	1;	%9 OCGT
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	9.967950000000002	2.0436	2.8	136	136	6135	10.0	2;	%9 biomass
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.654581100686952	0.0	0.0	350	350	3609	17.0	3;	%9 offwind
	10	0.0	0.0	0.0	-0.0	1.0	100.0	1	8.678973732622076	0.0	0.0	297	297	1615	38.0	4;	%9 onwind
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
	2	358618.2704780342	143447.3081912137	3	0.0	1829.6850534593582	0.0	0.0	486.09501821981223;	%0 CCGT
	2	483782.21921682166	193512.88768672865	3	0.0	2468.276628657253	0.0	0.0	221.35890241441408;	%0 OCGT
	2	3274.188301792864	1309.6753207171457	3	0.0	640.8667648840994	0.0	0.0	1365.0196130490774;	%0 biomass
	2	0.0	0.0	3	0.0	1.0624557994322381	0.0	0.0	944.5010173554292;	%0 offwind
	2	0.0	0.0	3	0.0	1.056821268851245	0.0	0.0	471.65447315736884;	%0 onwind
	2	0.0	0.0	3	0.0	0.8457051591370095	0.0	0.0	174.75902551182548;	%0 solar
	2	360229.2283891943	144091.6913556777	3	0.0	1837.904226475481	0.0	0.0	486.09501821981223;	%1 CCGT
	2	463327.22033845156	185330.8881353806	3	0.0	2363.9143894818953	0.0	0.0	221.35890241441408;	%1 OCGT
	2	3274.1193538885454	1309.6477415554182	3	0.0	640.8532695025534	0.0	0.0	1365.0196130490774;	%1 biomass
	2	0.0	0.0	3	0.0	1.075626761933853	0.0	0.0	471.65447315736884;	%1 onwind
	2	0.0	0.0	3	0.0	0.8656179647485729	0.0	0.0	174.75902551182548;	%1 solar
	2	343407.5051856011	137363.00207424044	3	0.0	1752.0791080898014	0.0	0.0	486.09501821981223;	%2 CCGT
	2	476063.3284395069	190425.33137580275	3	0.0	2428.8945328546265	0.0	0.0	221.35890241441408;	%2 OCGT
	2	3274.0923558232885	1309.6369423293154	3	0.0	640.8479850897022	0.0	0.0	1365.0196130490772;	%2 biomass
	2	0.0	0.0	3	0.0	1.0798596595619965	0.0	0.0	471.65447315736884;	%2 onwind
	2	0.0	0.0	3	0.0	0.8583143727650708	0.0	0.0	174.75902551182548;	%2 solar
	2	317980.22890705866	127192.09156282347	3	0.0	1622.3481066686666	0.0	0.0	486.09501821981223;	%3 CCGT
	2	456791.3843607563	182716.55374430254	3	0.0	2330.5682875548787	0.0	0.0	221.35890241441408;	%3 OCGT
	2	3274.1855701220584	1309.6742280488234	3	0.0	640.8662302059225	0.0	0.0	1365.0196130490774;	%3 biomass
	2	0.0	0.0	3	0.0	1.0812294105674807	0.0	0.0	967.0879585892485;	%3 offwind
	2	0.0	0.0	3	0.0	1.0618808132828672	0.0	0.0	471.6544731573689;	%3 onwind
	2	0.0	0.0	3	0.0	0.8535636267142581	0.0	0.0	174.75902551182548;	%3 solar
	2	340489.0043752764	136195.60175011054	3	0.0	1737.1887978330424	0.0	0.0	486.09501821981223;	%4 CCGT
	2	469927.06393932144	187970.82557572855	3	0.0	2397.587060914905	0.0	0.0	221.3589024144141;	%4 OCGT
	2	3274.1891185482814	1309.6756474193123	3	0.0	640.8669247501039	0.0	0.0	1365.0196130490774;	%4 biomass
	2	0.0	0.0	3	0.0	1.0647138940504528	0.0	0.0	471.65447315736884;	%4 onwind
	2	0.0	0.0	3	0.0	0.840882244684098	0.0	0.0	174.75902551182548;	%4 solar
	2	354510.5491792271	141804.21967169084	3	0.0	1808.7272917307503	0.0	0.0	486.09501821981223;	%5 CCGT
	2	472105.85306197515	188842.34122479009	3	0.0	2408.7033319488532	0.0	0.0	221.35890241441408;	%5 OCGT
	2	3274.204881737991	1309.6819526951965	3	0.0	640.8700101268332	0.0	0.0	1365.0196130490774;	%5 biomass
	2	0.0	0.0	3	0.0	1.075757421576654	0.0	0.0	471.65447315736884;	%5 onwind
	2	0.0	0.0	3	0.0	0.870387542101847	0.0	0.0	174.75902551182548;	%5 solar
	2	337187.7996181659	134875.11984726635	3	0.0	1720.3459164192136	0.0	0.0	486.09501821981223;	%6 CCGT
	2	476287.5856407154	190515.03425628616	3	0.0	2430.0387022485475	0.0	0.0	221.3589024144141;	%6 OCGT
	2	3274.159029846737	1309.6636119386947	3	0.0	640.861035397678	0.0	0.0	1365.0196130490774;	%6 biomass
	2	0.0	0.0	3	0.0	1.0908428858655497	0.0	0.0	471.65447315736884;	%6 onwind
	2	0.0	0.0	3	0.0	0.8686733249953484	0.0	0.0	174.75902551182546;	%6 solar
	2	342675.2691023908	137070.1076409563	3	0.0	1748.3432097060752	0.0	0.0	486.09501821981223;	%7 CCGT
	2	479465.51279324753	191786.20511729902	3	0.0	2446.252616292079	0.0	0.0	221.35890241441408;	%7 OCGT
	2	3274.2019430255423	1309.6807772102168	3	0.0	640.8694349237702	0.0	0.0	1365.0196130490774;	%7 biomass
	2	0.0	0.0	3	0.0	1.082403822651966	0.0	0.0	471.65447315736884;	%7 onwind
	2	0.0	0.0	3	0.0	0.8533955604998925	0.0	0.0	174.75902551182548;	%7 solar
	2	368543.0864232628	147417.23456930512	3	0.0	1880.321869506443	0.0	0.0	486.09501821981223;	%8 CCGT
	2	482102.66892572894	192841.06757029158	3	0.0	2459.7074945190247	0.0	0.0	221.35890241441408;	%8 OCGT
	2	3274.2274636812185	1309.6909854724872	3	0.0	640.8744301587822	0.0	0.0	1365.0196130490774;	%8 biomass
	2	0.0	0.0	3	0.0	1.0896908750655965	0.0	0.0	471.65447315736884;	%8 onwind
	2	0.0	0.0	3	0.0	0.857831538648998	0.0	0.0	174.75902551182548;	%8 solar
	2	371149.1409374151	148459.65637496606	3	0.0	1893.6180660072198	0.0	0.0	486.09501821981223;	%9 CCGT
	2	483917.1162634672	193566.84650538687	3	0.0	2468.9648788952404	0.0	0.0	221.3589024144141;	%9 OCGT
	2	3274.098605655283	1309.639442262113	3	0.0	640.8492083881939	0.0	0.0	1365.0196130490774;	%9 biomass
	2	0.0	0.0	3	0.0	1.0652952052244127	0.0	0.0	998.8034187377012;	%9 offwind
	2	0.0	0.0	3	0.0	1.0740683906844435	0.0	0.0	471.65447315736884;	%9 onwind
	2	0.0	0.0	3	0.0	0.8596297333255521	0.0	0.0	174.75902551182548;	%9 solar
	2	0.0	0.0	3	0.0	0.41572773784153316	0.0	0.0	7625844.314637915;	%DE1 0 PHS
	2	0.0	0.0	3	0.0	0.31179580338114987	0.0	0.0	7625844.314637915;	%DE1 0 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.45292744474548785	0.0	0.0	7625844.314637915;	%DE1 12 PHS
	2	0.0	0.0	3	0.0	0.3396955835591159	0.0	0.0	7625844.314637915;	%DE1 12 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4726055010764821	0.0	0.0	7625844.314637915;	%DE1 14 PHS
	2	0.0	0.0	3	0.0	0.3544541258073616	0.0	0.0	7625844.314637915;	%DE1 14 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3946730961078202	0.0	0.0	7625844.314637915;	%DE1 17 PHS
	2	0.0	0.0	3	0.0	0.29600482208086515	0.0	0.0	7625844.314637915;	%DE1 17 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4695465280618309	0.0	0.0	7625844.314637915;	%DE1 2 PHS
	2	0.0	0.0	3	0.0	0.35215989604637316	0.0	0.0	7625844.314637915;	%DE1 2 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.45170407781485156	0.0	0.0	7625844.314637915;	%DE1 20 PHS
	2	0.0	0.0	3	0.0	0.33877805836113867	0.0	0.0	7625844.314637915;	%DE1 20 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4299623226553939	0.0	0.0	7625844.314637915;	%DE1 21 PHS
	2	0.0	0.0	3	0.0	0.32247174199154544	0.0	0.0	7625844.314637915;	%DE1 21 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4364439792324717	0.0	0.0	7625844.314637915;	%DE1 25 hydro
	2	0.0	0.0	3	0.0	0.4315316564815954	0.0	0.0	7625844.314637915;	%DE1 27 PHS
	2	0.0	0.0	3	0.0	0.32364874236119656	0.0	0.0	7625844.314637915;	%DE1 27 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.42021471890451234	0.0	0.0	7625844.314637915;	%DE1 31 PHS
	2	0.0	0.0	3	0.0	0.31516103917838423	0.0	0.0	7625844.314637915;	%DE1 31 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4124910372543152	0.0	0.0	7625844.314637915;	%DE1 33 PHS
	2	0.0	0.0	3	0.0	0.3093682779407364	0.0	0.0	7625844.314637915;	%DE1 33 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.395020587489003	0.0	0.0	7625844.314637915;	%DE1 36 PHS
	2	0.0	0.0	3	0.0	0.29626544061675225	0.0	0.0	7625844.314637915;	%DE1 36 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4356954396199975	0.0	0.0	7625844.314637915;	%DE1 43 PHS
	2	0.0	0.0	3	0.0	0.32677157971499815	0.0	0.0	7625844.314637915;	%DE1 43 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3999471247192728	0.0	0.0	7625844.314637915;	%DE1 49 PHS
	2	0.0	0.0	3	0.0	0.2999603435394546	0.0	0.0	7625844.314637915;	%DE1 49 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4023532886686612	0.0	0.0	7625844.314637915;	%DE1 50 PHS
	2	0.0	0.0	3	0.0	0.3017649665014959	0.0	0.0	7625844.314637915;	%DE1 50 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.46627630880229576	0.0	0.0	7625844.314637915;	%DE1 52 PHS
	2	0.0	0.0	3	0.0	0.3497072316017218	0.0	0.0	7625844.314637915;	%DE1 52 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4268397709985287	0.0	0.0	7625844.314637915;	%DE1 54 PHS
	2	0.0	0.0	3	0.0	0.3201298282488965	0.0	0.0	7625844.314637915;	%DE1 54 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.4436177810247967	0.0	0.0	7625844.314637915;	%DE1 61 PHS
	2	0.0	0.0	3	0.0	0.3327133357685975	0.0	0.0	7625844.314637915;	%DE1 61 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.44973881782104297	0.0	0.0	7625844.314637915;	%DE1 7 PHS
	2	0.0	0.0	3	0.0	0.33730411336578225	0.0	0.0	7625844.314637915;	%DE1 7 PHS (pump mode)
	2	0.0	0.0	3	0.0	0.3993149122905933	0.0	0.0	7625844.314637915;	%DE1 71 hydro
	2	0.0	0.0	3	0.0	0.46479831926442544	0.0	0.0	7625844.314637915;	%DE1 76 PHS
	2	0.0	0.0	3	0.0	0.3485987394483191	0.0	0.0	7625844.314637915;	%DE1 76 PHS (pump mode)
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

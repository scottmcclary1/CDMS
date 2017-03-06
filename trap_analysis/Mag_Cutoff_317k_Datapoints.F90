    SUBROUTINE mag_cutoff_317k (window_length, frequency, mag_cutoff)
    
    IMPLICIT NONE
    
    REAL(4), INTENT(OUT) :: mag_cutoff                         ! FFT magnitude analysis cutoff value
    REAL(4), INTENT(IN) :: frequency
    REAL(4) :: a
    REAL(4) :: b
    REAL(4) :: c

    INTEGER, INTENT(IN) :: window_length 



    IF	(window_length	==	1000	) THEN
	    a	=	15646.99643	
	    b	=	22219.00388	
	    c	=	7913.81399	
    ELSE IF	(window_length	==	2000	) THEN
	    a	=	20777.07995	
	    b	=	29096.27929	
	    c	=	13570.87033	
    ELSE IF	(window_length	==	3000	) THEN
	    a	=	25576.67232	
	    b	=	28106.33457	
	    c	=	13354.62459	
    ELSE IF	(window_length	==	4000	) THEN
	    a	=	29342.76079	
	    b	=	28596.24703	
	    c	=	13810.87678	
    ELSE IF	(window_length	==	5000	) THEN
	    a	=	32635.73234	
	    b	=	29214.86217	
	    c	=	14305.33317	
    ELSE IF	(window_length	==	6000	) THEN
	    a	=	35843.10261	
	    b	=	28991.02106	
	    c	=	14202.28151	
    ELSE IF	(window_length	==	7000	) THEN
	    a	=	38790.77738	
	    b	=	28811.54061	
	    c	=	14117.88769	
    ELSE IF	(window_length	==	8000	) THEN
	    a	=	41487.63652	
	    b	=	28781.856	
	    c	=	14137.65115	
    ELSE IF	(window_length	==	9000	) THEN
	    a	=	43998.30838	
	    b	=	28780.71677	
	    c	=	14154.81802	
    ELSE IF	(window_length	==	10000	) THEN
	    a	=	46346.54556	
	    b	=	28811.95759	
	    c	=	14180.07718	
    ELSE IF	(window_length	==	11000	) THEN
	    a	=	48623.42972	
	    b	=	28737.61249	
	    c	=	14129.53684	
    ELSE IF	(window_length	==	12000	) THEN
	    a	=	50828.32039	
	    b	=	28623.09507	
	    c	=	14053.20007	
    ELSE IF	(window_length	==	13000	) THEN
	    a	=	52931.05691	
	    b	=	28558.0702	
	    c	=	14013.94036	
    ELSE IF	(window_length	==	14000	) THEN
	    a	=	54937.37287	
	    b	=	28531.36103	
	    c	=	13999.20225	
    ELSE IF	(window_length	==	15000	) THEN
	    a	=	56857.01628	
	    b	=	28526.69517	
	    c	=	13993.98761	
    ELSE IF	(window_length	==	16000	) THEN
	    a	=	58712.92635	
	    b	=	28523.22215	
	    c	=	13994.82524	
    ELSE IF	(window_length	==	17000	) THEN
	    a	=	60518.43144	
	    b	=	28503.95272	
	    c	=	13982.37933	
    ELSE IF	(window_length	==	18000	) THEN
	    a	=	62265.34085	
	    b	=	28504.94647	
	    c	=	13981.95385	
    ELSE IF	(window_length	==	19000	) THEN
	    a	=	63975.31013	
	    b	=	28491.59943	
	    c	=	13969.98391	
    ELSE IF	(window_length	==	20000	) THEN
	    a	=	65644.49545	
	    b	=	28478.03735	
	    c	=	13957.06547	
    ELSE IF	(window_length	==	21000	) THEN
	    a	=	67264.71686	
	    b	=	28485.74819	
	    c	=	13960.8601	
    ELSE IF	(window_length	==	22000	) THEN
	    a	=	68845.32274	
	    b	=	28498.80826	
	    c	=	13968.11338	
    ELSE IF	(window_length	==	23000	) THEN
	    a	=	70374.48952	
	    b	=	28540.75508	
	    c	=	13995.08146	
    ELSE IF	(window_length	==	24000	) THEN
	    a	=	71866.78842	
	    b	=	28587.4523	
	    c	=	14024.98342	
    ELSE IF	(window_length	==	25000	) THEN
	    a	=	73327.55118	
	    b	=	28633.44447	
	    c	=	14054.20055	
    ELSE IF	(window_length	==	26000	) THEN
	    a	=	74765.52216	
	    b	=	28666.01113	
	    c	=	14073.88742	
    ELSE IF	(window_length	==	27000	) THEN
	    a	=	76181.76409	
	    b	=	28688.95738	
	    c	=	14087.82707	
    ELSE IF	(window_length	==	28000	) THEN
	    a	=	77577.79933	
	    b	=	28700.31128	
	    c	=	14094.71601	
    ELSE IF	(window_length	==	29000	) THEN
	    a	=	78953.03908	
	    b	=	28704.20403	
	    c	=	14096.93309	
    ELSE IF	(window_length	==	30000	) THEN
	    a	=	80302.80756	
	    b	=	28709.87228	
	    c	=	14100.30626	
    ELSE IF	(window_length	==	31000	) THEN
	    a	=	81633.00871	
	    b	=	28711.25409	
	    c	=	14101.50922	
    ELSE IF	(window_length	==	32000	) THEN
	    a	=	82941.34969	
	    b	=	28712.30495	
	    c	=	14102.94973	
    ELSE IF	(window_length	==	33000	) THEN
	    a	=	84227.44241	
	    b	=	28714.6306	
	    c	=	14106.06194	
    ELSE IF	(window_length	==	34000	) THEN
	    a	=	85495.16854	
	    b	=	28716.64357	
	    c	=	14108.05499	
    ELSE IF	(window_length	==	35000	) THEN
	    a	=	86746.04498	
	    b	=	28718.34673	
	    c	=	14111.50542	
    ELSE IF	(window_length	==	36000	) THEN
	    a	=	87982.06371	
	    b	=	28715.39332	
	    c	=	14112.235	
    ELSE IF	(window_length	==	37000	) THEN
	    a	=	89202.01197	
	    b	=	28712.01994	
	    c	=	14112.76193	
    ELSE IF	(window_length	==	38000	) THEN
	    a	=	90407.49455	
	    b	=	28707.34972	
	    c	=	14112.34596	
    ELSE IF	(window_length	==	39000	) THEN
	    a	=	91591.52841	
	    b	=	28709.00935	
	    c	=	14115.95905	
    ELSE IF	(window_length	==	40000	) THEN
	    a	=	92759.51997	
	    b	=	28711.89533	
	    c	=	14120.89165	
    ELSE IF	(window_length	==	41000	) THEN
	    a	=	93914.04349	
	    b	=	28713.2713	
	    c	=	14124.18019	
    ELSE IF	(window_length	==	42000	) THEN
	    a	=	95054.29781	
	    b	=	28715.08294	
	    c	=	14128.41688	
    ELSE IF	(window_length	==	43000	) THEN
	    a	=	96181.01046	
	    b	=	28718.79076	
	    c	=	14134.09649	
    ELSE IF	(window_length	==	44000	) THEN
	    a	=	97293.17855	
	    b	=	28722.1216	
	    c	=	14138.76413	
    ELSE IF	(window_length	==	45000	) THEN
	    a	=	98392.74697	
	    b	=	28726.07964	
	    c	=	14144.0901	
    ELSE IF	(window_length	==	46000	) THEN
	    a	=	99480.77252	
	    b	=	28730.07443	
	    c	=	14150.72927	
    ELSE IF	(window_length	==	47000	) THEN
	    a	=	100556.0308	
	    b	=	28734.20512	
	    c	=	14156.55766	
    ELSE IF	(window_length	==	48000	) THEN
	    a	=	101616.8917	
	    b	=	28740.82995	
	    c	=	14163.22619	
    ELSE IF	(window_length	==	49000	) THEN
	    a	=	102667.4506	
	    b	=	28747.14266	
	    c	=	14169.76331	
    ELSE IF	(window_length	==	50000	) THEN
	    a	=	103707.889	
	    b	=	28751.2001	
	    c	=	14173.77372	
    ELSE IF	(window_length	==	51000	) THEN
	    a	=	104733.0873	
	    b	=	28762.49191	
	    c	=	14183.93759	
    ELSE IF	(window_length	==	52000	) THEN
	    a	=	105746.6941	
	    b	=	28774.71462	
	    c	=	14194.9036	
    ELSE IF	(window_length	==	53000	) THEN
	    a	=	106750.1302	
	    b	=	28787.34068	
	    c	=	14204.79853	
    ELSE IF	(window_length	==	54000	) THEN
	    a	=	107744.5523	
	    b	=	28798.63482	
	    c	=	14213.4372	
    ELSE IF	(window_length	==	55000	) THEN
	    a	=	108725.6371	
	    b	=	28814.84421	
	    c	=	14225.93467	
    ELSE IF	(window_length	==	56000	) THEN
	    a	=	109697.9462	
	    b	=	28830.51844	
	    c	=	14237.1442	
    ELSE IF	(window_length	==	57000	) THEN
	    a	=	110660.2632	
	    b	=	28847.34879	
	    c	=	14249.77162	
    ELSE IF	(window_length	==	58000	) THEN
	    a	=	111614.6527	
	    b	=	28863.81075	
	    c	=	14261.92683	
    ELSE IF	(window_length	==	59000	) THEN
	    a	=	112557.6667	
	    b	=	28882.24904	
	    c	=	14275.25553	
    ELSE IF	(window_length	==	60000	) THEN
	    a	=	113490.1322	
	    b	=	28901.95986	
	    c	=	14288.69861	
    ELSE IF	(window_length	==	61000	) THEN
	    a	=	114412.7503	
	    b	=	28923.42959	
	    c	=	14303.2503	
    ELSE IF	(window_length	==	62000	) THEN
	    a	=	115335.8635	
	    b	=	28935.55072	
	    c	=	14310.88005	
    ELSE IF	(window_length	==	63000	) THEN
	    a	=	116248.372	
	    b	=	28949.34709	
	    c	=	14319.66953	
    ELSE IF	(window_length	==	64000	) THEN
	    a	=	117158.0678	
	    b	=	28959.05678	
	    c	=	14326.19171	
    ELSE IF	(window_length	==	65000	) THEN
	    a	=	118060.4763	
	    b	=	28967.73349	
	    c	=	14331.11004	
    ELSE IF	(window_length	==	66000	) THEN
	    a	=	118955.7824	
	    b	=	28976.48805	
	    c	=	14336.5402	
    ELSE IF	(window_length	==	67000	) THEN
	    a	=	119843.5779	
	    b	=	28986.13646	
	    c	=	14342.37352	
    ELSE IF	(window_length	==	68000	) THEN
	    a	=	120723.3237	
	    b	=	28996.12711	
	    c	=	14348.47666	
    ELSE IF	(window_length	==	69000	) THEN
        a	=	121594.0774	
        b	=	29008.95121	
        c	=	14356.6472
    ELSE IF	(window_length	==	70000	) THEN
	    a	=	122464.0545	
	    b	=	29015.42163	
	    c	=	14360.21255	
    ELSE IF	(window_length	==	71000	) THEN
	    a	=	123325.8817	
	    b	=	29024.03838	
	    c	=	14365.73238	
    ELSE IF	(window_length	==	72000	) THEN
	    a	=	124181.2624	
	    b	=	29031.81012	
	    c	=	14370.42819	
    ELSE IF	(window_length	==	73000	) THEN
	    a	=	125034.5526	
	    b	=	29036.31171	
	    c	=	14373.23708	
    ELSE IF	(window_length	==	74000	) THEN
	    a	=	125879.363	
	    b	=	29042.69885	
	    c	=	14376.96159	
    ELSE IF	(window_length	==	75000	) THEN
	    a	=	126719.2942	
	    b	=	29048.22326	
	    c	=	14380.09199	
    ELSE IF	(window_length	==	76000	) THEN
	    a	=	127553.0331	
	    b	=	29054.76723	
	    c	=	14384.73697	
    ELSE IF	(window_length	==	77000	) THEN
	    a	=	128385.5823	
	    b	=	29056.96551	
	    c	=	14386.45921	
    ELSE IF	(window_length	==	78000	) THEN
	    a	=	129215.2898	
	    b	=	29056.67327	
	    c	=	14386.16097	
    ELSE IF	(window_length	==	79000	) THEN
	    a	=	130034.5211	
	    b	=	29062.1135	
	    c	=	14390.34135	
    ELSE IF	(window_length	==	80000	) THEN
	    a	=	130852.9857	
	    b	=	29063.38577	
	    c	=	14391.7683	
    ELSE IF	(window_length	==	81000	) THEN
	    a	=	131666.7223	
	    b	=	29063.982	
	    c	=	14392.72604	
    ELSE IF	(window_length	==	82000	) THEN
	    a	=	132473.5575	
	    b	=	29067.0706	
	    c	=	14395.47902	
    ELSE IF	(window_length	==	83000	) THEN
	    a	=	133276.1083	
	    b	=	29069.24344	
	    c	=	14397.537	
    ELSE IF	(window_length	==	84000	) THEN
	    a	=	134075.9616	
	    b	=	29069.47187	
	    c	=	14398.05586	
    ELSE IF	(window_length	==	85000	) THEN
	    a	=	134870.156	
	    b	=	29070.54192	
	    c	=	14399.01585	
    ELSE IF	(window_length	==	86000	) THEN
	    a	=	135657.4319	
	    b	=	29073.27049	
	    c	=	14401.46766	
    ELSE IF	(window_length	==	87000	) THEN
	    a	=	136440.8606	
	    b	=	29075.95279	
	    c	=	14404.06626	
    ELSE IF	(window_length	==	88000	) THEN
	    a	=	137221.9968	
	    b	=	29075.20064	
	    c	=	14403.76662	
    ELSE IF	(window_length	==	89000	) THEN
	    a	=	137998.1805	
	    b	=	29076.30493	
	    c	=	14405.33244	
    ELSE IF	(window_length	==	90000	) THEN
	    a	=	138768.3657	
	    b	=	29078.73065	
	    c	=	14407.65072	
    ELSE IF	(window_length	==	91000	) THEN
	    a	=	139532.8219	
	    b	=	29082.08361	
	    c	=	14410.80012	
    ELSE IF	(window_length	==	92000	) THEN
	    a	=	140293.3555	
	    b	=	29085.0541	
	    c	=	14413.52213	
    ELSE IF	(window_length	==	93000	) THEN
	    a	=	141050.2566	
	    b	=	29087.1578	
	    c	=	14415.47865	
    ELSE IF	(window_length	==	94000	) THEN
	    a	=	141803.7422	
	    b	=	29088.68849	
	    c	=	14417.16254	
    ELSE IF	(window_length	==	95000	) THEN
	    a	=	142553.434	
	    b	=	29089.39875	
	    c	=	14417.83303	
    ELSE IF	(window_length	==	96000	) THEN
	    a	=	143298.5982	
	    b	=	29090.82022	
	    c	=	14418.86539	
    ELSE IF	(window_length	==	97000	) THEN
	    a	=	144038.7586	
	    b	=	29092.89078	
	    c	=	14420.33055	
    ELSE IF	(window_length	==	98000	) THEN
	    a	=	144774.8209	
	    b	=	29095.8036	
	    c	=	14422.24809	
    ELSE IF	(window_length	==	99000	) THEN
	    a	=	145508.0499	
	    b	=	29097.15291	
	    c	=	14423.05513	
    ELSE IF	(window_length	==	100000	) THEN
	    a	=	146235.6355	
	    b	=	29100.36321	
	    c	=	14425.43149	
    ELSE IF	(window_length	==	101000	) THEN
	    a	=	146963.654	
	    b	=	29100.38642	
	    c	=	14425.79888	
    ELSE IF	(window_length	==	102000	) THEN
	    a	=	147686.0762	
	    b	=	29101.98505	
	    c	=	14426.74278	
    ELSE IF	(window_length	==	103000	) THEN
	    a	=	148408.1178	
	    b	=	29100.23569	
	    c	=	14425.0712	
    ELSE IF	(window_length	==	104000	) THEN
	    a	=	149123.7799	
	    b	=	29101.47899	
	    c	=	14425.66822	
    ELSE IF	(window_length	==	105000	) THEN
	    a	=	149836.1558	
	    b	=	29102.71447	
	    c	=	14426.10409	
    ELSE IF	(window_length	==	106000	) THEN
	    a	=	150546.7595	
	    b	=	29102.90945	
	    c	=	14426.08311	
    ELSE IF	(window_length	==	107000	) THEN
	    a	=	151253.4155	
	    b	=	29102.84775	
	    c	=	14425.5549	
    ELSE IF	(window_length	==	108000	) THEN
	    a	=	151957.566	
	    b	=	29102.99402	
	    c	=	14425.46689	
    ELSE IF	(window_length	==	109000	) THEN
	    a	=	152657.595	
	    b	=	29104.45107	
	    c	=	14426.34878	
    ELSE IF	(window_length	==	110000	) THEN
	    a	=	153353.4294	
	    b	=	29106.11629	
	    c	=	14427.34224	
    ELSE IF	(window_length	==	111000	) THEN
	    a	=	154048.6888	
	    b	=	29105.80935	
	    c	=	14426.56583	
    ELSE IF	(window_length	==	112000	) THEN
	    a	=	154740.7665	
	    b	=	29105.09048	
	    c	=	14425.6633	
    ELSE IF	(window_length	==	113000	) THEN
	    a	=	155431.0102	
	    b	=	29104.227	
	    c	=	14424.79153	
    ELSE IF	(window_length	==	114000	) THEN
	    a	=	156114.2713	
	    b	=	29106.28624	
	    c	=	14425.87621	
    ELSE IF	(window_length	==	115000	) THEN
	    a	=	156795.3707	
	    b	=	29108.0575	
	    c	=	14426.88177	
    ELSE IF	(window_length	==	116000	) THEN
	    a	=	157473.5664	
	    b	=	29110.01712	
	    c	=	14428.19474	
    ELSE IF	(window_length	==	117000	) THEN
	    a	=	158151.1396	
	    b	=	29109.96735	
	    c	=	14428.00456	
    ELSE IF	(window_length	==	118000	) THEN
	    a	=	158825.3752	
	    b	=	29110.25421	
	    c	=	14427.98481	
    ELSE IF	(window_length	==	119000	) THEN
	    a	=	159497.7042	
	    b	=	29109.58946	
	    c	=	14427.2693	
    ELSE IF	(window_length	==	120000	) THEN
	    a	=	160168.0691	
	    b	=	29108.8971	
	    c	=	14426.74158	
    ELSE IF	(window_length	==	121000	) THEN
	    a	=	160835.9082	
	    b	=	29107.96022	
	    c	=	14425.79565	
    ELSE IF	(window_length	==	122000	) THEN
	    a	=	161500.5736	
	    b	=	29107.07051	
	    c	=	14425.05202	
    ELSE IF	(window_length	==	123000	) THEN
	    a	=	162164.4993	
	    b	=	29104.59942	
	    c	=	14423.0269	
    ELSE IF	(window_length	==	124000	) THEN
	    a	=	162826.4437	
	    b	=	29102.26452	
	    c	=	14421.24226	
    ELSE IF	(window_length	==	125000	) THEN
	    a	=	163484.3094	
	    b	=	29100.84482	
	    c	=	14420.07523	
    ELSE IF	(window_length	==	126000	) THEN
	    a	=	164137.7092	
	    b	=	29101.31021	
	    c	=	14420.49556	
    ELSE IF	(window_length	==	127000	) THEN
	    a	=	164788.6316	
	    b	=	29100.88181	
	    c	=	14420.24406	
    ELSE IF	(window_length	==	128000	) THEN
	    a	=	165438.4914	
	    b	=	29099.41911	
	    c	=	14419.03968	
    ELSE IF	(window_length	==	129000	) THEN
	    a	=	166085.1319	
	    b	=	29098.674	
	    c	=	14418.6969	
    ELSE IF	(window_length	==	130000	) THEN
	    a	=	166728.646	
	    b	=	29098.56965	
	    c	=	14418.5703	
    ELSE IF	(window_length	==	131000	) THEN
	    a	=	167367.7523	
	    b	=	29100.30833	
	    c	=	14420.08542	
    ELSE IF	(window_length	==	132000	) THEN
	    a	=	168009.4255	
	    b	=	29098.19085	
	    c	=	14418.89143	
    ELSE IF	(window_length	==	133000	) THEN
	    a	=	168646.188	
	    b	=	29097.58963	
	    c	=	14418.69429	
    ELSE IF	(window_length	==	134000	) THEN
	    a	=	169279.4095	
	    b	=	29098.28149	
	    c	=	14419.30683	
    ELSE IF	(window_length	==	135000	) THEN
	    a	=	169911.1823	
	    b	=	29097.67112	
	    c	=	14419.06657	
    ELSE IF	(window_length	==	136000	) THEN
	    a	=	170542.2242	
	    b	=	29096.04455	
	    c	=	14418.17084	
    ELSE IF	(window_length	==	137000	) THEN
	    a	=	171169.6373	
	    b	=	29095.92429	
	    c	=	14418.61004	
    ELSE IF	(window_length	==	138000	) THEN
	    a	=	171793.9376	
	    b	=	29096.51577	
	    c	=	14419.485	
    ELSE IF	(window_length	==	139000	) THEN
	    a	=	172417.9435	
	    b	=	29095.18971	
	    c	=	14418.66202	
    ELSE IF	(window_length	==	140000	) THEN
	    a	=	173039.185	
	    b	=	29094.37169	
	    c	=	14418.27948	
    ELSE IF	(window_length	==	141000	) THEN
	    a	=	173659.3058	
	    b	=	29092.69579	
	    c	=	14417.45986	
    ELSE IF	(window_length	==	142000	) THEN
	    a	=	174276.4071	
	    b	=	29091.99799	
	    c	=	14417.42901	
    ELSE IF	(window_length	==	143000	) THEN
	    a	=	174891.065	
	    b	=	29090.86817	
	    c	=	14416.82089	
    ELSE IF	(window_length	==	144000	) THEN
	    a	=	175502.4512	
	    b	=	29091.13054	
	    c	=	14417.30369	
    ELSE IF	(window_length	==	145000	) THEN
	    a	=	176113.4067	
	    b	=	29089.99567	
	    c	=	14416.81601	
    ELSE IF	(window_length	==	146000	) THEN
	    a	=	176721.4794	
	    b	=	29089.38241	
	    c	=	14416.63662	
    ELSE IF	(window_length	==	147000	) THEN
	    a	=	177329.7086	
	    b	=	29087.52303	
	    c	=	14415.72578	
    ELSE IF	(window_length	==	148000	) THEN
	    a	=	177934.0236	
	    b	=	29087.04354	
	    c	=	14415.87939	
    ELSE IF	(window_length	==	149000	) THEN
	    a	=	178534.8997	
	    b	=	29087.37382	
	    c	=	14416.43719	
    ELSE IF	(window_length	==	150000	) THEN
	    a	=	179134.7945	
	    b	=	29087.07885	
	    c	=	14416.59506	
    ELSE IF	(window_length	==	151000	) THEN
	    a	=	179732.2547	
	    b	=	29086.87643	
	    c	=	14416.76931	
    ELSE IF	(window_length	==	152000	) THEN
	    a	=	180329.1246	
	    b	=	29085.8049	
	    c	=	14416.39171	
    ELSE IF	(window_length	==	153000	) THEN
	    a	=	180922.4455	
	    b	=	29085.60234	
	    c	=	14416.33953	
    ELSE IF	(window_length	==	154000	) THEN
	    a	=	181514.5397	
	    b	=	29085.41741	
	    c	=	14416.67869	
    ELSE IF	(window_length	==	155000	) THEN
	    a	=	182104.6939	
	    b	=	29085.26034	
	    c	=	14416.85047	
    ELSE IF	(window_length	==	156000	) THEN
	    a	=	182691.32	
	    b	=	29085.93246	
	    c	=	14417.69549	
    ELSE IF	(window_length	==	157000	) THEN
	    a	=	183277.319	
	    b	=	29085.53438	
	    c	=	14417.80282	
    ELSE IF	(window_length	==	158000	) THEN
	    a	=	183860.2821	
	    b	=	29086.14263	
	    c	=	14418.4835	
    ELSE IF	(window_length	==	159000	) THEN
	    a	=	184443.2509	
	    b	=	29085.78297	
	    c	=	14418.53016	
    ELSE IF	(window_length	==	160000	) THEN
	    a	=	185024.0203	
	    b	=	29085.56004	
	    c	=	14418.69349	
    ELSE IF	(window_length	==	161000	) THEN
	    a	=	185601.8784	
	    b	=	29085.68012	
	    c	=	14419.09476	
    ELSE IF	(window_length	==	162000	) THEN
	    a	=	186178.992	
	    b	=	29085.42388	
	    c	=	14419.3168	
    ELSE IF	(window_length	==	163000	) THEN
	    a	=	186752.8262	
	    b	=	29086.02503	
	    c	=	14420.26659	
    ELSE IF	(window_length	==	164000	) THEN
	    a	=	187327.1969	
	    b	=	29084.97058	
	    c	=	14419.85546	
    ELSE IF	(window_length	==	165000	) THEN
	    a	=	187898.8159	
	    b	=	29084.73663	
	    c	=	14420.03336	
    ELSE IF	(window_length	==	166000	) THEN
	    a	=	188467.138	
	    b	=	29085.40802	
	    c	=	14420.92079	
    ELSE IF	(window_length	==	167000	) THEN
	    a	=	189035.5203	
	    b	=	29084.91747	
	    c	=	14420.91988	
    ELSE IF	(window_length	==	168000	) THEN
	    a	=	189601.9453	
	    b	=	29084.58031	
	    c	=	14420.90521	
    ELSE IF	(window_length	==	169000	) THEN
	    a	=	190167.05	
	    b	=	29083.72869	
	    c	=	14420.55446	
    ELSE IF	(window_length	==	170000	) THEN
	    a	=	190730.3874	
	    b	=	29083.09724	
	    c	=	14420.317	
    ELSE IF	(window_length	==	171000	) THEN
	    a	=	191292.0153	
	    b	=	29082.46814	
	    c	=	14420.04786	
    ELSE IF	(window_length	==	172000	) THEN
	    a	=	191852.8634	
	    b	=	29081.61	
	    c	=	14419.76753	
    ELSE IF	(window_length	==	173000	) THEN
	    a	=	192410.888	
	    b	=	29081.15931	
	    c	=	14419.58776	
    ELSE IF	(window_length	==	174000	) THEN
	    a	=	192967.548	
	    b	=	29080.78931	
	    c	=	14419.58147	
    ELSE IF	(window_length	==	175000	) THEN
	    a	=	193523.4544	
	    b	=	29079.72483	
	    c	=	14418.94522	
    ELSE IF	(window_length	==	176000	) THEN
	    a	=	194076.6697	
	    b	=	29079.21617	
	    c	=	14418.75885	
    ELSE IF	(window_length	==	177000	) THEN
	    a	=	194629.747	
	    b	=	29078.12332	
	    c	=	14418.19053	
    ELSE IF	(window_length	==	178000	) THEN
	    a	=	195181.764	
	    b	=	29076.45486	
	    c	=	14417.19111	
    ELSE IF	(window_length	==	179000	) THEN
	    a	=	195729.4048	
	    b	=	29076.70769	
	    c	=	14417.42382	
    ELSE IF	(window_length	==	180000	) THEN
	    a	=	196278.1604	
	    b	=	29075.27221	
	    c	=	14416.52328	
    ELSE IF	(window_length	==	181000	) THEN
	    a	=	196824.1867	
	    b	=	29074.50199	
	    c	=	14416.10505	
    ELSE IF	(window_length	==	182000	) THEN
	    a	=	197368.6719	
	    b	=	29073.3048	
	    c	=	14415.34113	
    ELSE IF	(window_length	==	183000	) THEN
	    a	=	197912.3365	
	    b	=	29072.42187	
	    c	=	14414.93	
    ELSE IF	(window_length	==	184000	) THEN
	    a	=	198452.1877	
	    b	=	29072.30459	
	    c	=	14414.82065	
    ELSE IF	(window_length	==	185000	) THEN
	    a	=	198993.5076	
	    b	=	29070.71438	
	    c	=	14413.77475	
    ELSE IF	(window_length	==	186000	) THEN
	    a	=	199531.0757	
	    b	=	29070.49003	
	    c	=	14413.75262	
    ELSE IF	(window_length	==	187000	) THEN
	    a	=	200068.8024	
	    b	=	29068.80574	
	    c	=	14412.49718	
    ELSE IF	(window_length	==	188000	) THEN
	    a	=	200603.8954	
	    b	=	29067.89343	
	    c	=	14411.84701	
    ELSE IF	(window_length	==	189000	) THEN
	    a	=	201139.4656	
	    b	=	29066.29073	
	    c	=	14410.78646	
    ELSE IF	(window_length	==	190000	) THEN
	    a	=	201672.3815	
	    b	=	29065.17689	
	    c	=	14409.95516	
    ELSE IF	(window_length	==	191000	) THEN
	    a	=	202205.1332	
	    b	=	29063.41822	
	    c	=	14408.75078	
    ELSE IF	(window_length	==	192000	) THEN
	    a	=	202734.8396	
	    b	=	29062.36476	
	    c	=	14407.8722	
    ELSE IF	(window_length	==	193000	) THEN
	    a	=	203263.9425	
	    b	=	29061.06632	
	    c	=	14406.92214	
    ELSE IF	(window_length	==	194000	) THEN
	    a	=	203791.0403	
	    b	=	29060.07706	
	    c	=	14406.16071	
    ELSE IF	(window_length	==	195000	) THEN
	    a	=	204316.3341	
	    b	=	29059.72993	
	    c	=	14405.99696	
    ELSE IF	(window_length	==	196000	) THEN
	    a	=	204840.3336	
	    b	=	29059.08793	
	    c	=	14405.54995	
    ELSE IF	(window_length	==	197000	) THEN
	    a	=	205363.9086	
	    b	=	29057.55551	
	    c	=	14404.27402	
    ELSE IF	(window_length	==	198000	) THEN
	    a	=	205885.6825	
	    b	=	29057.02033	
	    c	=	14404.14476	
    ELSE IF	(window_length	==	199000	) THEN
	    a	=	206405.2136	
	    b	=	29056.40739	
	    c	=	14403.57366	
    ELSE IF	(window_length	==	200000	) THEN
	    a	=	206922.3461	
	    b	=	29056.62987	
	    c	=	14403.62375	
    ELSE IF	(window_length	==	201000	) THEN
	    a	=	207438.316	
	    b	=	29056.74393	
	    c	=	14403.48471	
    ELSE IF	(window_length	==	202000	) THEN
	    a	=	207954.3135	
	    b	=	29056.13384	
	    c	=	14402.98203	
    ELSE IF	(window_length	==	203000	) THEN
	    a	=	208469.0386	
	    b	=	29055.55188	
	    c	=	14402.50754	
    ELSE IF	(window_length	==	204000	) THEN
	    a	=	208982.4972	
	    b	=	29055.09122	
	    c	=	14402.15452	
    ELSE IF	(window_length	==	205000	) THEN
	    a	=	209492.9489	
	    b	=	29055.51484	
	    c	=	14402.46953	
    ELSE IF	(window_length	==	206000	) THEN
	    a	=	210002.8891	
	    b	=	29055.66403	
	    c	=	14402.45326	
    ELSE IF	(window_length	==	207000	) THEN
	    a	=	210511.7514	
	    b	=	29055.58629	
	    c	=	14402.27833	
    ELSE IF	(window_length	==	208000	) THEN
	    a	=	211020.2942	
	    b	=	29054.89853	
	    c	=	14401.73532	
    ELSE IF	(window_length	==	209000	) THEN
	    a	=	211526.2992	
	    b	=	29054.88545	
	    c	=	14401.63924	
    ELSE IF	(window_length	==	210000	) THEN
	    a	=	212029.384	
	    b	=	29056.08407	
	    c	=	14402.38437	
    ELSE IF	(window_length	==	211000	) THEN
	    a	=	212531.6149	
	    b	=	29057.07504	
	    c	=	14403.05171	
    ELSE IF	(window_length	==	212000	) THEN
	    a	=	213032.6427	
	    b	=	29058.06024	
	    c	=	14403.66018	
    ELSE IF	(window_length	==	213000	) THEN
	    a	=	213534.3381	
	    b	=	29057.89785	
	    c	=	14403.51993	
    ELSE IF	(window_length	==	214000	) THEN
	    a	=	214034.5798	
	    b	=	29057.87028	
	    c	=	14403.40673	
    ELSE IF	(window_length	==	215000	) THEN
	    a	=	214532.1587	
	    b	=	29058.83019	
	    c	=	14404.01853	
    ELSE IF	(window_length	==	216000	) THEN
	    a	=	215028.0403	
	    b	=	29059.70108	
	    c	=	14404.37412	
    ELSE IF	(window_length	==	217000	) THEN
	    a	=	215522.089	
	    b	=	29061.22548	
	    c	=	14405.24395	
    ELSE IF	(window_length	==	218000	) THEN
	    a	=	216018.4509	
	    b	=	29061.04739	
	    c	=	14405.07519	
    ELSE IF	(window_length	==	219000	) THEN
	    a	=	216511.4378	
	    b	=	29061.97118	
	    c	=	14405.64355	
    ELSE IF	(window_length	==	220000	) THEN
	    a	=	217003.8613	
	    b	=	29062.34357	
	    c	=	14405.73139	
    ELSE IF	(window_length	==	221000	) THEN
	    a	=	217494.7394	
	    b	=	29063.0284	
	    c	=	14406.00879	
    ELSE IF	(window_length	==	222000	) THEN
	    a	=	217984.8895	
	    b	=	29063.38814	
	    c	=	14406.04757	
    ELSE IF	(window_length	==	223000	) THEN
	    a	=	218473.5027	
	    b	=	29064.06917	
	    c	=	14406.33312	
    ELSE IF	(window_length	==	224000	) THEN
	    a	=	218961.1887	
	    b	=	29064.68864	
	    c	=	14406.61385	
    ELSE IF	(window_length	==	225000	) THEN
	    a	=	219447.2572	
	    b	=	29065.66848	
	    c	=	14407.07326	
    ELSE IF	(window_length	==	226000	) THEN
	    a	=	219931.8884	
	    b	=	29066.87169	
	    c	=	14407.74475	
    ELSE IF	(window_length	==	227000	) THEN
	    a	=	220418.0143	
	    b	=	29066.72392	
	    c	=	14407.5597	
    ELSE IF	(window_length	==	228000	) THEN
	    a	=	220901.3059	
	    b	=	29067.61905	
	    c	=	14408.17966	
    ELSE IF	(window_length	==	229000	) THEN
	    a	=	221382.903	
	    b	=	29068.94623	
	    c	=	14409.06613	
    ELSE IF	(window_length	==	230000	) THEN
	    a	=	221863.3807	
	    b	=	29070.13591	
	    c	=	14409.82157	
    ELSE IF	(window_length	==	231000	) THEN
	    a	=	222344.2874	
	    b	=	29070.30415	
	    c	=	14409.72889	
    ELSE IF	(window_length	==	232000	) THEN
	    a	=	222823.9965	
	    b	=	29070.63259	
	    c	=	14409.78986	
    ELSE IF	(window_length	==	233000	) THEN
	    a	=	223301.4994	
	    b	=	29071.64138	
	    c	=	14410.42441	
    ELSE IF	(window_length	==	234000	) THEN
	    a	=	223778.1544	
	    b	=	29072.59979	
	    c	=	14410.99695	
    ELSE IF	(window_length	==	235000	) THEN
	    a	=	224252.8677	
	    b	=	29073.98169	
	    c	=	14411.77887	
    ELSE IF	(window_length	==	236000	) THEN
	    a	=	224728.2294	
	    b	=	29074.62563	
	    c	=	14412.13225	
    ELSE IF	(window_length	==	237000	) THEN
	    a	=	225200.5093	
	    b	=	29076.06354	
	    c	=	14412.99325	
    ELSE IF	(window_length	==	238000	) THEN
	    a	=	225674.1367	
	    b	=	29076.4633	
	    c	=	14413.18472	
    ELSE IF	(window_length	==	239000	) THEN
	    a	=	226147.1892	
	    b	=	29076.59762	
	    c	=	14413.20596	
    ELSE IF	(window_length	==	240000	) THEN
	    a	=	226617.0562	
	    b	=	29077.94979	
	    c	=	14414.04715	
    ELSE IF	(window_length	==	241000	) THEN
	    a	=	227088.6847	
	    b	=	29077.80322	
	    c	=	14413.80223	
    ELSE IF	(window_length	==	242000	) THEN
	    a	=	227557.9425	
	    b	=	29078.49259	
	    c	=	14414.19128	
    ELSE IF	(window_length	==	243000	) THEN
	    a	=	228026.8169	
	    b	=	29079.08094	
	    c	=	14414.66665	
    ELSE IF	(window_length	==	244000	) THEN
	    a	=	228494.881	
	    b	=	29079.18838	
	    c	=	14414.48483	
    ELSE IF	(window_length	==	245000	) THEN
	    a	=	228962.0107	
	    b	=	29079.24773	
	    c	=	14414.33361	
    ELSE IF	(window_length	==	246000	) THEN
	    a	=	229427.0919	
	    b	=	29079.98652	
	    c	=	14414.75909	
    ELSE IF	(window_length	==	247000	) THEN
	    a	=	229892.2364	
	    b	=	29080.1307	
	    c	=	14414.77371	
    ELSE IF	(window_length	==	248000	) THEN
	    a	=	230356.3827	
	    b	=	29080.19224	
	    c	=	14414.77281	
    ELSE IF	(window_length	==	249000	) THEN
	    a	=	230818.4371	
	    b	=	29081.12018	
	    c	=	14415.28919	
    ELSE IF	(window_length	==	250000	) THEN
	    a	=	231281.9121	
	    b	=	29080.67043	
	    c	=	14414.79778	
    ELSE IF	(window_length	==	251000	) THEN
	    a	=	231741.6308	
	    b	=	29081.89772	
	    c	=	14415.60113	
    ELSE IF	(window_length	==	252000	) THEN
	    a	=	232201.7056	
	    b	=	29082.27055	
	    c	=	14415.81968	
    ELSE IF	(window_length	==	253000	) THEN
	    a	=	232663.0664	
	    b	=	29081.61267	
	    c	=	14415.22561	
    ELSE IF	(window_length	==	254000	) THEN
	    a	=	233122.6006	
	    b	=	29081.23354	
	    c	=	14414.8315	
    ELSE IF	(window_length	==	255000	) THEN
	    a	=	233579.7225	
	    b	=	29081.92537	
	    c	=	14415.18368	
    ELSE IF	(window_length	==	256000	) THEN
	    a	=	234036.0535	
	    b	=	29082.38713	
	    c	=	14415.36652	
    ELSE IF	(window_length	==	257000	) THEN
	    a	=	234491.0949	
	    b	=	29083.06211	
	    c	=	14415.83334	
    ELSE IF	(window_length	==	258000	) THEN
	    a	=	234945.7134	
	    b	=	29083.51063	
	    c	=	14416.15665	
    ELSE IF	(window_length	==	259000	) THEN
	    a	=	235398.8107	
	    b	=	29084.37854	
	    c	=	14416.67302	
    ELSE IF	(window_length	==	260000	) THEN
	    a	=	235849.273	
	    b	=	29086.10684	
	    c	=	14417.80942	
    ELSE IF	(window_length	==	261000	) THEN
	    a	=	236301.5132	
	    b	=	29086.30536	
	    c	=	14417.85798	
    ELSE IF	(window_length	==	262000	) THEN
	    a	=	236752.8389	
	    b	=	29086.36926	
	    c	=	14417.82396	
    ELSE IF	(window_length	==	263000	) THEN
	    a	=	237201.5004	
	    b	=	29087.60757	
	    c	=	14418.65049	
    ELSE IF	(window_length	==	264000	) THEN
	    a	=	237650.227	
	    b	=	29088.27642	
	    c	=	14418.96353	
    ELSE IF	(window_length	==	265000	) THEN
	    a	=	238097.8494	
	    b	=	29089.2292	
	    c	=	14419.48276	
    ELSE IF	(window_length	==	266000	) THEN
	    a	=	238545.3622	
	    b	=	29089.93426	
	    c	=	14420.01079	
    ELSE IF	(window_length	==	267000	) THEN
	    a	=	238991.1637	
	    b	=	29091.00449	
	    c	=	14420.71806	
    ELSE IF	(window_length	==	268000	) THEN
	    a	=	239436.3939	
	    b	=	29091.95506	
	    c	=	14421.37688	
    ELSE IF	(window_length	==	269000	) THEN
	    a	=	239881.7876	
	    b	=	29092.40787	
	    c	=	14421.63371	
    ELSE IF	(window_length	==	270000	) THEN
	    a	=	240324.1836	
	    b	=	29093.81857	
	    c	=	14422.55274	
    ELSE IF	(window_length	==	271000	) THEN
	    a	=	240767.6381	
	    b	=	29094.14818	
	    c	=	14422.72825	
    ELSE IF	(window_length	==	272000	) THEN
	    a	=	241208.2638	
	    b	=	29095.53503	
	    c	=	14423.71048	
    ELSE IF	(window_length	==	273000	) THEN
	    a	=	241650.7527	
	    b	=	29095.85387	
	    c	=	14424.09091	
    ELSE IF	(window_length	==	274000	) THEN
	    a	=	242090.3288	
	    b	=	29096.99776	
	    c	=	14424.80878	
    ELSE IF	(window_length	==	275000	) THEN
	    a	=	242531.6259	
	    b	=	29096.8265	
	    c	=	14424.45796	
    ELSE IF	(window_length	==	276000	) THEN
	    a	=	242970.4183	
	    b	=	29097.47859	
	    c	=	14424.74585	
    ELSE IF	(window_length	==	277000	) THEN
	    a	=	243407.5074	
	    b	=	29098.80605	
	    c	=	14425.80501	
    ELSE IF	(window_length	==	278000	) THEN
	    a	=	243845.026	
	    b	=	29099.42261	
	    c	=	14426.27991	
    ELSE IF	(window_length	==	279000	) THEN
	    a	=	244281.0576	
	    b	=	29100.19688	
	    c	=	14426.68349	
    ELSE IF	(window_length	==	280000	) THEN
	    a	=	244716.2385	
	    b	=	29101.01097	
	    c	=	14427.06675	
    ELSE IF	(window_length	==	281000	) THEN
	    a	=	245150.5606	
	    b	=	29102.02824	
	    c	=	14427.68385	
    ELSE IF	(window_length	==	282000	) THEN
	    a	=	245584.5034	
	    b	=	29102.76262	
	    c	=	14428.32112	
    ELSE IF	(window_length	==	283000	) THEN
	    a	=	246017.1844	
	    b	=	29103.84616	
	    c	=	14429.13123	
    ELSE IF	(window_length	==	284000	) THEN
	    a	=	246450.0919	
	    b	=	29104.44773	
	    c	=	14429.52762	
    ELSE IF	(window_length	==	285000	) THEN
	    a	=	246880.7787	
	    b	=	29105.64841	
	    c	=	14430.2684	
    ELSE IF	(window_length	==	286000	) THEN
	    a	=	247312.2204	
	    b	=	29106.13768	
	    c	=	14430.51912	
    ELSE IF	(window_length	==	287000	) THEN
	    a	=	247742.2966	
	    b	=	29107.01626	
	    c	=	14431.30867	
    ELSE IF	(window_length	==	288000	) THEN
	    a	=	248172.4423	
	    b	=	29107.68054	
	    c	=	14431.93018	
    ELSE IF	(window_length	==	289000	) THEN
	    a	=	248601.0485	
	    b	=	29108.21847	
	    c	=	14432.15811	
    ELSE IF	(window_length	==	290000	) THEN
	    a	=	249028.7981	
	    b	=	29109.2048	
	    c	=	14432.74375	
    ELSE IF	(window_length	==	291000	) THEN
	    a	=	249454.436	
	    b	=	29110.72363	
	    c	=	14433.63994	
    ELSE IF	(window_length	==	292000	) THEN
	    a	=	249881.1105	
	    b	=	29111.45619	
	    c	=	14434.33292	
    ELSE IF	(window_length	==	293000	) THEN
	    a	=	250306.6027	
	    b	=	29112.23627	
	    c	=	14435.0069	
    ELSE IF	(window_length	==	294000	) THEN
	    a	=	250731.4292	
	    b	=	29113.19132	
	    c	=	14435.66071	
    ELSE IF	(window_length	==	295000	) THEN
	    a	=	251155.689	
	    b	=	29114.09399	
	    c	=	14436.07852	
    ELSE IF	(window_length	==	296000	) THEN
	    a	=	251581.484	
	    b	=	29113.70241	
	    c	=	14435.64008	
    ELSE IF	(window_length	==	297000	) THEN
	    a	=	252001.6546	
	    b	=	29115.86007	
	    c	=	14437.46729	
    ELSE IF	(window_length	==	298000	) THEN
	    a	=	252422.5538	
	    b	=	29117.04457	
	    c	=	14438.49578	
    ELSE IF	(window_length	==	299000	) THEN
	    a	=	252845.0493	
	    b	=	29117.40596	
	    c	=	14438.80296	
    ELSE IF	(window_length	==	300000	) THEN
	    a	=	253265.8691	
	    b	=	29118.2553	
	    c	=	14439.09999	
    ELSE IF	(window_length	==	301000	) THEN
	    a	=	253686.897	
	    b	=	29118.52218	
	    c	=	14438.88852	
    ELSE IF	(window_length	==	302000	) THEN
	    a	=	254104.7457	
	    b	=	29120.10578	
	    c	=	14440.3167	
    ELSE IF	(window_length	==	303000	) THEN
	    a	=	254521.2848	
	    b	=	29121.80204	
	    c	=	14441.88798	
    ELSE IF	(window_length	==	304000	) THEN
	    a	=	254940.6655	
	    b	=	29122.05675	
	    c	=	14442.0925	
    ELSE IF	(window_length	==	305000	) THEN
	    a	=	255361.3993	
	    b	=	29121.07776	
	    c	=	14440.6903	
    ELSE IF	(window_length	==	306000	) THEN
	    a	=	255776.2625	
	    b	=	29122.81614	
	    c	=	14441.5466	
    ELSE IF	(window_length	==	307000	) THEN
	    a	=	256192.1094	
	    b	=	29123.43919	
	    c	=	14442.33633	
    ELSE IF	(window_length	==	308000	) THEN
	    a	=	256606.0102	
	    b	=	29124.79065	
	    c	=	14444.17264	
    ELSE IF	(window_length	==	309000	) THEN
	    a	=	257020.6361	
	    b	=	29125.51584	
	    c	=	14444.36465	
    ELSE IF	(window_length	==	310000	) THEN
	    a	=	257441.6288	
	    b	=	29122.96004	
	    c	=	14440.61161	
    ELSE IF	(window_length	==	311000	) THEN
	    a	=	257854.5288	
	    b	=	29124.09472	
	    c	=	14439.92987	
    ELSE IF	(window_length	==	312000	) THEN
	    a	=	258271.4682	
	    b	=	29122.3104	
	    c	=	14439.39176	
    ELSE IF	(window_length	==	313000	) THEN
	    a	=	258674.0537	
	    b	=	29127.25205	
	    c	=	14445.5134	
    ELSE IF	(window_length	==	314000	) THEN
	    a	=	259090.2714	
	    b	=	29125.44746	
	    c	=	14442.05926	
    ELSE IF	(window_length	==	315000	) THEN
	    a	=	259528.8304	
	    b	=	29113.94902	
	    c	=	14419.68477	
    ELSE IF	(window_length	==	316000	) THEN
	    a	=	260005.7794	
	    b	=	29084.31259	
	    c	=	14367.99183	
    ELSE IF	(window_length	==	317000	) THEN
	    a	=	260594.6185	
	    b	=	29006.81551	
	    c	=	14229.86991	
    END IF
          
    
    mag_cutoff = (a * EXP(b / (frequency + c))) + 30000
    
    
    END SUBROUTINE  mag_cutoff_317k
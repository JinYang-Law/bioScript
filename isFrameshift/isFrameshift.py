import bisect

import pandas
import pandas as pd
from collections import namedtuple

'''
input
ETV6:NM_001987.4:exon4-NTRK3:NM_001012338.2:exon12

refgene
84	NM_001987.4	chr12	+	11802787	12048325	11803061	12043980	8	11802787,11905383,11992073,12006360,12022357,12037378,12038859,12043874,	11803094,11905513,11992238,12006495,12022903,12037521,12038960,12048325,	0	ETV6	cmpl	cmpl	0,0,1,1,1,1,0,2,
27	NM_152263.3	chr1	-	154134288	154164611	154140412	154164494	10	154134288,154141780,154142875,154143124,154143888,154145383,154145559,154148590,154163661,154164377,	154140416,154141859,154142945,154143187,154143964,154145454,154145677,154148724,154163787,154164611,	0	TPM3	cmpl	cmpl	2,1,0,0,2,0,2,0,0,0,
831	NM_004521.3	chr10	-	32297942	32345353	32304456	32344901	26	32297942,32304436,32306070,32306979,32307243,32307429,32308785,32309949,32310153,32311067,32311775,32317355,32320000,32321633,32322772,32323617,32324449,32324817,32326181,32326447,32327090,32327705,32328254,32329311,32337391,32344775,	32300444,32304587,32306287,32307084,32307315,32307490,32308887,32310059,32310215,32311185,32311964,32317499,32320207,32321702,32322966,32323766,32324595,32324922,32326306,32326535,32327146,32327754,32328359,32329385,32337479,32345353,	0	KIF5B	cmpl	cmpl	-1,1,0,0,0,2,2,0,1,0,0,0,0,0,1,2,0,0,1,0,1,0,0,1,0,0,
157	NM_001012338.2	chr15	-	88419987	88799962	88420165	88799384	20	88419987,88423500,88428924,88472421,88476242,88483853,88576087,88669501,88670392,88671941,88678331,88679129,88679697,88680634,88690565,88726648,88727455,88799136,88799514,88799874,	88420351,88423659,88428966,88472665,88476415,88483984,88576276,88669604,88670457,88671965,88678628,88679271,88679840,88680792,88690634,88726720,88727530,88799399,88799717,88799962,	0	NTRK3	cmpl	cmpl	0,0,0,2,0,1,1,0,1,1,1,0,1,2,2,2,2,0,-1,-1,
1781	NM_002529.3	chr1	+	156830670	156851642	156830726	156851434	17	156830670,156834145,156834519,156836701,156837895,156838296,156841414,156843424,156844174,156844362,156844697,156845311,156845871,156846191,156848913,156849790,156851248,	156830938,156834220,156834591,156836770,156838041,156838439,156841547,156843751,156844192,156844418,156844800,156845458,156846002,156846364,156849154,156849949,156851642,	0	NTRK1	cmpl	cmpl	0,2,2,2,2,1,0,1,1,1,0,1,1,0,2,0,0,

output
  GeneF       TransF   ExonF  GeneD  ... ExonnD up_cdsLen down_cdsLen mod
0  ETV6  NM_001987.4   exon4  NTRK3  ...     12       463        1292   0
1  TPM3  NM_152263.3  exon10  NTRK1  ...      9      6982        1214   0

'''


class Frameshift(object):
    def __init__(self, ref, fusion, output):
        self.ref = ref
        self.fusion = fusion
        self.out = output

    def _read_ref(self):
        ref_Gene = pd.read_csv(self.ref, sep="\t",
                               names=['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
                                      'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat',
                                      'cdsEndStat', 'exonFrames'])
        return ref_Gene


    def _read_fusion(self):
        fus = pd.read_csv(self.fusion, sep="-|:", names=["GeneF", 'TransF', 'ExonF', "GeneD", 'TransD', 'ExonD'])
        fus['ExonnF'] = fus['ExonF'].map(lambda x: x.replace("exon", ""))
        fus['ExonnD'] = fus['ExonD'].map(lambda x: x.replace("exon", ""))
        return fus


    def _upstream_parser(self, series, detail):
        starts = list(filter(None, detail['exonStarts'].split(",")))
        ends = list(filter(None, detail['exonEnds'].split(",")))
        Exonlen = [[starts[ind], ends[ind]] for ind in range(len(starts))]
        start_keys = [int(tp[0]) for tp in Exonlen]
        end_keys = [int(tp[1]) for tp in Exonlen]
        cds = 0

        if detail['strand'] == "+":
            ind = bisect.bisect_right(start_keys, int(detail['cdsStart']))
            if ind == 0:
                real_Len = Exonlen[ind:series['ExonnF']]
                cds = sum(map(lambda x: x[1] - x[0], real_Len))
            else:
                Exonlen[ind - 1][0] = detail['cdsStart']
                real_Len = Exonlen[ind - 1:int(series['ExonnF'])]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
        else:
            ind = bisect.bisect_left(end_keys, int(detail['cdsEnd']))
            if ind == len(end_keys):
                real_Len = Exonlen[series['ExonnF']:ind]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
            else:
                Exonlen[ind][1] = detail['cdsEnd']
                real_Len = Exonlen[len(end_keys) - int(series['ExonnF']):ind + 1]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))

        return cds


    def _downstream_parser(self, series, detail):
        starts = list(filter(None, detail['exonStarts'].split(",")))
        ends = list(filter(None, detail['exonEnds'].split(",")))
        Exonlen = [[starts[ind], ends[ind]] for ind in range(len(starts))]
        start_keys = [int(tp[0]) for tp in Exonlen]
        end_keys = [int(tp[1]) for tp in Exonlen]
        cds = 0

        if detail['strand'] == "+":
            ind = bisect.bisect_left(end_keys, int(detail['cdsEnd']))
            if ind == len(end_keys):
                real_Len = Exonlen[int(series['ExonnD'])-1:ind]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
            else:
                Exonlen[ind][1] = detail['cdsEnd']
                real_Len = Exonlen[int(series['ExonnD'])-1:ind+1]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
        else:
            ind = bisect.bisect_right(start_keys, int(detail['cdsStart']))
            if ind == 0:
                real_Len = Exonlen[ind:int(detail['exonCount']) - int(series['ExonnD']) + 1]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
            else:
                Exonlen[ind - 1][0] = detail['cdsStart']
                real_Len = Exonlen[ind - 1: int(detail['exonCount']) - int(series['ExonnD']) + 1]
                cds = sum(map(lambda x: int(x[1]) - int(x[0]), real_Len))
            pass

        return cds


    def _isShift(self, series, ref):
        ## 返回series


        ### 对上游融合基因进行分析
        up_info = ref[(ref['name'] == series['TransF'])]
        if up_info.empty: return False
        up_cdsLen = self._upstream_parser(series, up_info.iloc[0,])
        up_mod = up_cdsLen % 3

        ### 对下游融合基因进行分析
        down_info = ref[(ref['name'] == series['TransD'])]
        if down_info.empty: return False
        down_cdsLen = self._downstream_parser(series, down_info.iloc[0,])
        down_mod = down_cdsLen % 3
        mod = (up_mod + down_mod) % 3

        return (pandas.Series([up_cdsLen, down_cdsLen, mod], index=['up_cdsLen', 'down_cdsLen', 'mod']))


    def isFrame_shift(self):
        ### 读取refGene 文件

        ref = self._read_ref()

        ### 读取融合文件
        fus = self._read_fusion()

        ### 应用apply 函数， 遍历每一行进行处理
        ## 待优化，如果没有返回值判断
        shift_flags = fus.apply(self._isShift, axis=1, ref=ref)

        return (pd.concat([fus,shift_flags],axis=1))



def main():
    refGene = 'refGene'
    fusionList = 'fusionList'
    out = 'out'
    frameshift = Frameshift(refGene, fusionList, out)
    fshift_out = frameshift.isFrame_shift()
    print(fshift_out)

if __name__ == '__main__':
    main()

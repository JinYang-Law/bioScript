import bisect

import pandas
import pandas as pd
from collections import namedtuple


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

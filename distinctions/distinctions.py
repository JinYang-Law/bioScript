import os
import pandas as pd
from collections import namedtuple


class Distinguish(object):
    def __init__(self, panel, dir, suffix, out_prefix):
        self.panel = panel
        self.dir = dir
        self.suffix = suffix
        self.prefix = out_prefix
        pass

    @property
    def _walk(self):
        Path = namedtuple('Path', ['dirname', 'fname'])
        for (pth, dirs, fnames) in os.walk(self.dir):
            for fname in fnames:
                if fname.endswith(self.suffix):
                    yield (Path(pth, fname))  ## 对象 + 属性

    def _evaluate(self, series, bed):
        ## 对数据框的每一行进行判断,判断标准： 1. 染色体一样， 变异位点在bed 区域内
        ## 对每个Series执行结果后，会将结果整合在一起返回（若想有返回值，定义函数时需要return相应的值）
        flags = bed[(bed['Chr'] == series['Chr']) & (bed['Start'] <= series['Start']) & (bed['End'] >= series['Start'])]
        if flags.empty:
            return False
        else:
            return True

    def _single_validate(self, file, bed) -> object:
        variation = pd.read_csv(file, sep="\t", header=0)
        if variation.empty:
            return variation

        variation_Sed = variation.loc[:, ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Freq', 'Depth_US']]
        variation_Sed['Flags'] = variation_Sed.apply(self._evaluate, axis=1, bed=bed)
        return variation_Sed

    def _variation_write(self, info):
        info["Chr_Pos"] = info['Chr'] + "_" + info['Start'].astype("str") + "_" + \
                          info['End'].astype("str") + "_" + info['Ref'] + "_" + \
                          info['Alt']

        info_sorted = info.sort_values("Chr_Pos")
        info_sorted['nPos'] = info_sorted.groupby("Chr_Pos")['Chr_Pos'].transform('count')
        info_sorted.to_csv(self.prefix + "_pos.sorted.txt", sep="\t", index=False)

        info_sorted_bySample = info.sort_values("Sample")
        info_sorted_bySample['Sample_nPos'] = info_sorted_bySample.groupby("Sample")['Sample'].transform('count')
        info_sorted_bySample.to_csv(self.prefix + "_sample.sorted.txt", sep="\t", index=False)

        info_sorted_redup = info_sorted.drop_duplicates(subset=['Chr_Pos'], keep='first')
        info_sorted_redup = info_sorted_redup.loc[:, ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Freq', 'nPos']]
        info_sorted_redup.to_csv(self.prefix + "_pos.redup.txt", sep="\t", index=False)

    def get_files(self):
        for Path in self._walk:  # 循环对象，取出属性
            isfile = os.path.join(Path.dirname, Path.fname)
            if os.path.exists(isfile):
                print(isfile)

    def distinguish(self):
        ## 读取panel 文件; (待优化：表头，列数的判断)
        bed = pd.read_csv(self.panel, sep="\t", header=None, names=['Chr', 'Start', 'End'])

        ## 读取变异检测文件，判断变异位点是否发生在panel的范围内
        ## (待优化：变异检测文件是否存在判断）
        ## filter() 函数 None 的判断
        ## 返回数据框不要返回None

        variation_pos = []
        variation_neg = []
        for Path in self._walk:
            file = os.path.join(Path.dirname, Path.fname)
            variation_valid = self._single_validate(file, bed)

            if variation_valid.empty:
                continue

            if variation_valid['Flags'].any():
                # variation_valid['Sample'] = Path.fname.split(".")[0]
                variation_valid.insert(0, 'Sample', Path.fname.split(".")[0])
                variation_pos.append(variation_valid[variation_valid['Flags']])
            else:
                variation_valid.insert(0, 'Sample', Path.fname.split(".")[0])
                variation_neg.append(variation_valid)

        if variation_pos:
            variation_pos_out = pd.concat(variation_pos)
            self._variation_write(variation_pos_out)

        if variation_neg:
            variation_neg_out = pd.concat(variation_neg)
            variation_neg_out.to_csv(self.prefix + "_neg.txt", sep="\t")


def main():
    ## 待优化：args
    bed = "LYN5.REF.bed"
    path = "20211012-Pair_Pos"
    suffix = "ann"
    prefix = 'Case_2014'
    lyn5gene = Distinguish(bed, path, suffix, prefix)
    lyn5gene.distinguish()


if __name__ == '__main__':
    main()

'''
定义类
创建对象
获取属性


# 创建一个User对象
user = User(name='kongxx', sex='male', age=21)
 
# 也可以通过一个list来创建一个User对象，这里注意需要使用"_make"方法
user = User._make(['kongxx', 'male', 21])
 
print user
# User(name='user1', sex='male', age=21)
 
# 获取用户的属性
print user.name
print user.sex
print user.age

'''


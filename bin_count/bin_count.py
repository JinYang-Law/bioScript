import bisect
import pandas as pd
from collections import defaultdict



class BinCount(object):
    """
    :param  list for bin
    :param  list for depths
    :return pandas.DataFrame
    """

    def __init__(self, bin, depths):
        self.bin = bin
        self.depths = depths
        pass


    def bin_count(self):
        """ 
            简单说明 count 定义字典的方式，count = {} , count = defaultdict(int) 差别
            1. count = {} 如果在修改字典时，key 不存在字典中的话，会报错：ValuesError; 
                改进： 可用if 判断key 是否存在，不存在的初始化
            2. 调用  defaultdict(int)； 直接默认初始值
        """
        count = defaultdict(int)
        dpct = {}
        for item in self.depths:
            ind = bisect.bisect_right(self.bin, item)      ## 返回右侧index， [1, 101, 201, 301, 401]  target=101 返回index = 2   （]  左边包含，右边不包含
            count[self.bin[ind-1]] += 1

        count = sorted(count.items(), key=lambda x:x[0])
        dp = [dp for dp, ct in count]
        ct = [ct for dp, ct in count]

        for ind in range(len(count)):
            if ind == len(count)-1:
                dpct[str(dp[ind]) + 'X-'] = ct[ind]
                continue
            dpct[str(dp[ind]) + 'X-' + str(dp[ind+1]-1) + 'X'] = ct[ind]
        return pd.DataFrame([dpct])


    def bin_accum(self):
        count = defaultdict(int)
        dpct = {}
        for item in self.depths:
            ind = bisect.bisect_right(self.bin, item)      ## 返回左侧index， [1, 5, 8, 11, 15]  target=5 返回index = 2
            count[self.bin[ind-1]] += 1

        count = sorted(count.items(), key=lambda x:x[0])
        dp = [dp for dp, ct in count]
        ct = [ct for dp, ct in count]

        for ind in range(len(count)):
            dpct[">=" + str(dp[ind])] = sum(ct[ind:])

        return pd.DataFrame([dict(count)])


def main():
    ## 待优化： args

    bin = [i for i in range(1, 1000, 100)]
    depths = [275, 431, 401, 401, 372, 424, 529, 497, 202, 167, 681, 203, 205, 735, 553, 348, 416, 486, 535, 535, 535, 273, 590, 447, 567, 588, 825, 151, 728, 243, 400, 415, 453, 629, 305, 628, 236, 382, 205, 264, 930, 313, 573, 234, 266, 287, 631, 313, 360, 544, 290, 530, 290, 530, 477, 303, 288, 778, 167, 611, 161, 181, 277, 669, 623, 840, 193, 187, 158, 663, 424, 524, 567, 462, 82, 299, 331, 438, 452, 418, 549, 243, 209, 553, 714, 610, 831, 395, 542, 720, 336, 685, 141, 440, 761, 334, 354, 542, 265, 412, 322, 404, 507, 414, 363, 443, 663, 709, 501, 306, 521, 381, 754, 307, 305, 326, 541, 413, 413, 217, 217]

    binC = BinCount(bin, depths)
    dp_count = binC.bin_count()


'''
  ## 进阶,咋么处理
    bin = [("1X-100X", 1), ('101X-200X',100) ... ]
'''


if __name__ == '__main__':
    main()




















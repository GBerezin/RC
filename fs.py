import numpy as np
import pandas as pd
import os
import time
import fclasses
import sclasses


class Start:
    def __init__(self, prop):
        self.prpr = prop.Properties()

    def strt(self, lst2, stg2, gb3, m, dat_file, load):
        start = time.time()
        os.system('cls')
        self.prpr.gb3 = gb3
        if lst2 != 0:
            self.prpr.lim_st = '2_st'
        else:
            self.prpr.lim_st = '1_st'
        self.prpr.prp(dat_file)
        self.prpr.chgprop()
        if m == 'frame':
            print('Program =RC_frame.v1.0= designed by G. Berezin in 2020')
            clss = fclasses
            u_stg1 = np.zeros(3)
        else:
            print('Program =RC_slab.v1.0= designed by G. Berezin in 2020')
            clss = sclasses
            u_stg1 = np.zeros(6)
        stfs = clss.Stiffness(self.prpr)
        lds = clss.Loads(load, self.prpr)
        lds.ld()
        clc = clss.Calc(stfs, self.prpr, u_stg1)
        itr = clss.Iteration(stfs, clc, self.prpr, self.prpr.dat)
        n = len(lds.F)
        eps = np.zeros((n, 2))
        for iF in range(0, n):
            Fg = lds.F[iF]
            itr.itrn(Fg, iF)
            if m == 'frame':
                eps[iF, 0] = round(min(itr.eps), 6)
                eps[iF, 1] = round(max(itr.eps), 6)
            else:
                min1 = min(itr.eps1)
                min2 = min(itr.eps2)
                min3 = min(itr.strain)
                max1 = max(itr.eps1)
                max2 = max(itr.eps2)
                max3 = max(itr.strain)
                eps[iF, 0] = round(min(min1, min2, min3), 6)
                eps[iF, 1] = round(max(max1, max2, max3), 6)
        result = pd.DataFrame(eps, columns=['eps_min', 'eps_max'])
        print(result)
        print('eps_min=', min(eps[:, 0]), 'L/C', np.where(
            result['eps_min'] == min(eps[:, 0]))[0])
        print('eps_max=', max(eps[:, 1]), 'L/C', np.where(
            result['eps_max'] == max(eps[:, 1]))[0])
        jt = round((time.time() - start), 4)
        print("Job time:", jt, "seconds")
        return stfs, self.prpr, result

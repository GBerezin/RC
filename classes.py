import crack
import os


class First:
    def __init__(self, prpr, lds, cnvr, itr, itr_b, clc, rsl, stfs, iF, dat_file, fname):
        self.prpr = prpr
        self.lds = lds
        self.cnvr = cnvr
        self.itr = itr
        self.itr_b = itr_b
        self.clc = clc
        self.rsl = rsl
        self.stfs = stfs
        self.iF = iF
        self.dat_file = dat_file
        self.fname = fname
        self.repname = fname+'.txt'
        self.resname = fname+'.csv'
        self.lds.ld()
        self.Fg = self.lds.F[self.iF]
        self.crc = crack.Crack(self.cnvr, self.dat_file)

    def calc(self, report):
        self.prpr.chgprop()
        print(self.prpr.lim_st, 'analysis:', file=report)
        self.itr.itrn(self.Fg, self.iF)
        print('Solution:', self.clc.mes, file=report)
        if self.clc.mes == 'ok':
            print('Results after', self.itr.i, 'iterations:', file=report)
            self.rsl.rslt1(self.itr)
            if self.itr.resmes != 'Convergence did not reached !':
                self.rsl.result1.to_csv(os.path.join(
                    self.prpr.dpath, self.resname), sep=';', index=False)
                print(self.rsl.result1, file=report)
            print(self.itr.resmes, file=report)
            self.cnvr.converg(self.Fg, self.itr)
            print(self.cnvr.cvrg, file=report)


class Second(First):
    def calc(self, report):
        self.prpr.lim_st = '2_st'
        self.prpr.prp(self.dat_file)
        self.lds.ld()
        self.crc.scrc(self.itr, self.itr_b, self.iF,
                      self.lds, self.prpr, self.rsl, report)
        print(self.crc.mes1, file=report)
        print(self.crc.mes2, file=report)
        print(self.crc.mes3, file=report)
        print('Solution:', self.clc.mes, file=report)
        print(self.rsl.result1, file=report)
        self.rsl.result1.to_csv(os.path.join(
            self.prpr.dpath, 'res2.csv'), sep=';', index=False)
        print(self.itr.resmes, file=report)
        print(self.cnvr.cvrg, file=report)
        if self.crc.Sscrc != 0.0:
            self.prpr.dur = 'long'
            self.lds.ld()
            self.crc.acr(self.lds, self.itr, self.clc, self.prpr)
            acrc1 = self.crc.acrc
            self.prpr.dur = 'short'
            self.crc.fi1 = 1.0
            self.crc.n = 2
            self.lds.ld()
            self.crc.acr(self.lds, self.itr, self.clc, self.prpr)
            acrc2 = self.crc.acrc
            self.prpr.dur = 'long'
            self.crc.fi1 = 1.0
            self.crc.n = 3
            self.lds.ld()
            self.crc.acr(self.lds, self.itr, self.clc, self.prpr)
            acrc3 = self.crc.acrc
            acrc = acrc1+acrc2-acrc3
            print('acrc=', round(acrc*1000, 3), 'mm', file=report)

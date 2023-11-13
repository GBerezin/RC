import numpy as np
import pandas as pd


class Crack:
    def __init__(self, cnvr, dat_file):
        self.fi1 = 1.4
        self.fi2 = 0.5
        self.fi3 = 1.0
        self.ncrc = 100
        self.cnvr = cnvr
        self.dat_file = dat_file
        self.n = 1
        self.rescrack = np.zeros((3, 13))

    def scrc(self, itr, itr_b, iF, lds, prpr, rsl, report):
        self.iF = iF
        self.report = report
        print(prpr.lim_st, 'analysis:', file=self.report)
        rng = np.linspace(lds.F[self.iF] / self.ncrc, lds.F[iF], self.ncrc)
        prpr.prp(self.dat_file)
        i = 1
        epsbt2 = prpr.epsbt2
        F_b = rng[0]
        for F in rng:
            itr.itrn(F, self.iF)
            self.cnvr.converg(F, itr)
            if max(itr.eps) > epsbt2:
                itr_b.itrn(F_b, self.iF)
                self.cnvr.converg(F_b, itr_b)
                self.ds = itr_b.d_s()
                self.Fcrc = itr_b.fcrc(F_b)
                self.As = np.sum(prpr.pr[(itr_b.Sig > 0) & (
                    prpr.pr['T'] == 'concrete')]['A'])
                self.Abt = np.sum(prpr.pr[(itr_b.Sig > 0) & (
                    prpr.pr['T'] == 'rebar')]['A'])
                self.mes1 = 'eps,crc=' + \
                    str(round(max(itr_b.eps), 7)) + \
                    ' at ' + str(i-1)+' % of L/C values'
                prpr.chgprop()
                print('Sig,crc calculation:', file=self.report)
                itr.itrn(F, self.iF)
                self.Sscrc = max(itr.Sig)
                self.mes2 = 'Sig,crc='+str(round(self.Sscrc, 3)) + \
                    ' MPa at '+str(i) + ' % of L/C values'
                self.mes3 = 'Results after '+str(itr.i) + ' iterations:'
                self.k = i-1
                break
            else:
                self.Fcrc = itr_b.fcrc(F_b)
                self.Sscrc = 0.0
                self.mes1 = 'No crack!'
                self.mes2 = 'Results after '+str(itr.i) + ' iterations:'
                self.mes3 = ''
                self.k = i
            i = i+1
            F_b = F
        prpr.prp(self.dat_file)
        itr_b.itrn(F_b, self.iF)
        rsl.rslt1(itr_b)
        itr.itrn(lds.F[iF], self.iF)
        self.Ss = max(lds.F[iF])

    def acr(self, lds, itr, clc, prpr):
        print('Crack width calculation:', file=self.report)
        if self.Sscrc != 0.0:
            print('Crack width', 'Phi1=', self.fi1,
                  prpr.dur, '-term loads', 'calculation:', file=self.report)
            F = lds.F[self.iF]
            itr.itrn(F, self.iF)
            print('Solution:', clc.mes, file=self.report)
            self.Ss = max(itr.Sig)
            self.Es = prpr.E[self.Ss == itr.Sig][0]
            n_F = len(F)
            if n_F == 3:
                if F[2] > 0.0:
                    self.fi3 = 1.2
            else:
                if F[0] > 0.0 or F[1] > 0.0:
                    self.fi3 = 1.2
            if clc.mes == 'ok':
                self.crci()
                self.rslt2(self)
                blankIndex=[''] * len(self.result2)
                self.result2.index=blankIndex
                print(self.result2, file=self.report)
                self.rescrack[self.n-1, :] = self.res2[0, :]

    def crci(self):
        if self.Sscrc <= self.Ss:
            self.psi_s = 1-0.8*self.Sscrc/self.Ss
        else:
            self.psi_s = 0.2
        l = 0.5 * self.Abt / self.As * self.ds
        low = max(10 * self.ds, 0.1)
        up = min(40 * self.ds, 0.4)
        if l < low:
            self.ls = low
        elif l > up:
            self.ls = up
        else:
            self.ls = l
        self.acrc = self.fi1*self.fi2*self.fi3*self.psi_s*self.Ss/self.Es*self.ls

    def rslt2(self, crc):
        Abt = round(crc.Abt*10000, 2)
        As = round(crc.As*10000, 2)
        ds = round(crc.ds*1000, 0)
        ls = round(crc.ls, 3)
        Phi1 = crc.fi1
        Phi2 = crc.fi2
        Phi3 = crc.fi3
        Psi_s = round(crc.psi_s, 3)
        Sscrc = round(crc.Sscrc, 3)
        Sig_s = round(crc.Ss, 3)
        Es = round(crc.Es, 1)
        acrc = round(crc.acrc*1000, 3)
        self.n = crc.n
        self.res2 = np.array([str(self.n), Abt, As, ds, ls, Phi1, Phi2, Phi3,
                              Psi_s, Sscrc, Sig_s, Es, acrc]).reshape(1, 13)
        self.result2 = pd.DataFrame(self.res2, columns=['crc', 'Abt,cm^2', 'As,cm^2', 'ds,mm', 'ls,m', 'Phi1',
                                                        'Phi2', 'Phi3', 'Psi,s', 'Sscrc,MPa', 'Ss,MPa', 'Es,MPa', 'acrc,mm'])

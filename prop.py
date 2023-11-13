import numpy as np
import pandas as pd
import os
dpath = os.path.join('.', 'data')
tpath = os.path.join('.', 'test')


class Properties:
    lim_st = '1_st'
    dur = 'short'
    gb3 = 1.0
    epsbt2 = 0.00015

    @staticmethod
    def sigma(e, e_2, e_0, e_1, e1, e0, e2, S_1, S1, Rc, Rt, E):
        if e_0 >= e >= e_2:
            S = Rc
        elif e_0 < e < e_1:
            S = ((1 - S_1 / Rc) * (e - e_1)/(e_0 - e_1) + S_1/Rc)*Rc
        elif e_1 <= e < 0.0:
            S = E * e
        elif 0.0 < e <= e1 and S1 != 0.0:
            S = E * e
        elif e1 < e < e0 and Rt != 0.0:
            S = ((1 - S1 / Rt) * (e - e1)/(e0 - e1) + S1 / Rt) * Rt
        elif e0 <= e <= e2:
            S = Rt
        else:
            S = 0.0
        return S

    def prp(self, dat_file):
        self.dpath = dpath
        self.mat = pd.read_csv('materials.csv', sep=';')
        self.dat = pd.read_csv(dat_file, sep=';')
        if self.lim_st == '1_st':
            gb = 1.0
            gbt = 1.0
            self.gs = 1.0
            self.gb3 = self.gb3
            if self.dur == 'long':
                gb1 = 0.9
            else:
                gb1 = 1.0
        else:
            gb1 = 1.0
            gb = 1.3
            gbt = 1.5
            self.gs = 1.15
            self.gb3 = 1.0
        self.grade = self.dat['Grade'].values
        self.pr = pd.merge(self.dat, self.mat)
        self.concr = self.pr[self.pr['T'] == 'concrete']
        self.steel = self.pr[self.pr['T'] == 'rebar']
        self.grade_b = self.concr['Grade'].values
        self.grade_s = self.steel['Grade'].values
        self.e_2 = self.pr['eps2'].values
        self.e_0 = self.pr['eps0'].values
        self.e_1 = self.pr['eps1'].values
        self.e1 = self.pr['epst1'].values
        self.e0 = self.pr['epst0'].values
        self.e2 = self.pr['epst2'].values
        self.Rc = pd.concat(
            [self.concr['Rc']*gb*gb1*self.gb3, self.steel['Rc']]).values
        self.S_1 = self.pr['Sc1'].values
        self.S1 = self.pr['St1'].values
        self.Rt = pd.concat(
            [self.concr['Rt']*gbt*gb1, self.steel['Rt'] * self.gs]).values
        self.E = self.pr['E'].values
        self.ce_2 = self.e_2[self.pr['T'] == 'concrete']
        self.ce_0 = self.e_0[self.pr['T'] == 'concrete']
        self.ce_1 = self.e_1[self.pr['T'] == 'concrete']
        self.ce1 = self.e1[self.pr['T'] == 'concrete']
        self.ce0 = self.e0[self.pr['T'] == 'concrete']
        self.ce2 = self.e2[self.pr['T'] == 'concrete']
        self.cS_1 = self.S_1[self.pr['T'] == 'concrete']
        self.cS1 = self.S1[self.pr['T'] == 'concrete']
        self.cRc = self.Rc[self.pr['T'] == 'concrete']
        self.cRt = self.Rt[self.pr['T'] == 'concrete']
        self.cE = self.E[self.pr['T'] == 'concrete']
        self.re_2 = self.e_2[self.pr['T'] == 'rebar']
        self.re_0 = self.e_0[self.pr['T'] == 'rebar']
        self.re_1 = self.e_1[self.pr['T'] == 'rebar']
        self.re1 = self.e1[self.pr['T'] == 'rebar']
        self.re0 = self.e0[self.pr['T'] == 'rebar']
        self.re2 = self.e2[self.pr['T'] == 'rebar']
        self.rS_1 = self.S_1[self.pr['T'] == 'rebar']
        self.rS1 = self.S1[self.pr['T'] == 'rebar']
        self.rRc = self.Rc[self.pr['T'] == 'rebar']
        self.rRt = self.Rt[self.pr['T'] == 'rebar']
        self.rE = self.E[self.pr['T'] == 'rebar']

    def chgprop(self):
        R0 = np.zeros(len(self.concr))
        S1_s = self.steel['St1'].values
        Rt_s = self.steel['Rt'].values * self.gs
        self.S1 = np.append(R0, S1_s)
        self.Rt = np.append(R0, Rt_s)
        self.cS1 = R0
        self.cRt = R0


class Sigma:
    @staticmethod
    def sigma(e, e_2, e_0, e_1, e1, e0, e2, S_1, S1, Rc, Rt, E):
        if e_0 >= e >= e_2:
            S = Rc
        elif e_0 < e < e_1:
            S = ((1 - S_1 / Rc) * (e - e_1)/(e_0 - e_1) + S_1/Rc)*Rc
        elif e_1 <= e < 0.0:
            S = E * e
        elif 0.0 < e <= e1 and S1 != 0.0:
            S = E * e
        elif e1 < e < e0 and Rt != 0.0:
            S = ((1 - S1 / Rt) * (e - e1)/(e0 - e1) + S1 / Rt) * Rt
        elif e0 <= e <= e2:
            S = Rt
        else:
            S = 0.0
        return S

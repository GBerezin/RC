import numpy as np
import pandas as pd
import pickle
import math
import os


class Iteration:
    def __init__(self, stfs, clc, prpr, dat):
        self.stfs = stfs
        self.clc = clc
        self.prpr = prpr
        self.dat = dat
        self.vV = np.vectorize(V.vi)

    def fcrc(self, F):
        Fcrc = pd.DataFrame(F.reshape(1, 3), columns=['Mx', 'My', 'N'])
        return Fcrc

    def d_s(self):
        ds = round(
            math.sqrt(4*max(self.dat['A'].values[max(self.Sig) == self.Sig])/math.pi), 3)
        return ds

    def itrn(self, Fg, iF):
        acc = 0.0000001
        self.stfs.d = np.ones(self.stfs.n)
        self.u = self.clc.calc(self.stfs.d, Fg)
        v = self.vV(self.clc.S, self.prpr.E, self.clc.e)
        i = 0
        du = 0.0001
        self.resmes = 'Convergence is reached !'
        while du >= acc:
            i = i+1
            if min(self.clc.e) < -0.0035 or max(self.clc.e) > 0.025:
                self.clc.e = np.zeros(self.stfs.n)
                self.resmes = 'Convergence L/C ' +  str(iF) + ' did not reached!!!'
                print(self.resmes)
                break
            self.stfs.d = v
            u_f = self.clc.calc(self.stfs.d, Fg)
            v = self.vV(self.clc.S, self.prpr.E, self.clc.e)
            du = np.max(abs(self.u - u_f))
            self.u = u_f
        self.eps = self.clc.e
        self.Sig = self.clc.S
        self.i = i


class V:
    @staticmethod
    def vi(S, E, e):
        if e != 0.0:
            v = S/E/e
        else:
            v = 1.0
        return v


class Stiffness:
    def __init__(self, prpr):
        self.prpr = prpr
        self.x = self.prpr.dat['Zx'].values
        self.y = self.prpr.dat['Zy'].values
        self.A = self.prpr.dat['A'].values
        self.n = len(self.x)
        self.Z = np.transpose(np.array([self.x, self.y, np.ones(self.n)]))
        self.v = np.ones(self.n)
        self.vsigma = np.vectorize(self.prpr.sigma)

    @property
    def d(self):
        E = self.prpr.E
        D = np.array([
            [sum(self.A*self.x**2*E*self.v), sum(self.A*self.x *
                                                 self.y*E*self.v), sum(self.A*self.x*E*self.v)],
            [sum(self.A*self.x*self.y*E*self.v), sum(self.A *
                                                     self.y**2*E*self.v), sum(self.A*self.y*E*self.v)],
            [sum(self.A*self.x*E*self.v),   sum(self.A*self.y*E*self.v),
             sum(self.A*E*self.v)]
        ])
        return D

    @d.setter
    def d(self, v):
        self.v = v


class Calc:
    def __init__(self, stfs, prpr, u_stg1):
        self.stfs = stfs
        self.prpr = prpr
        self.u_stg1 = u_stg1

    def calc(self, D, F):
        try:
            u = np.linalg.solve(D, F)
            self.mes = 'ok'
        except np.linalg.LinAlgError as var1:
            self.mes = var1
            D = np.eye(3)
            Fi = np.array(([0.0, 0.0, 0.0]))
            u = np.linalg.solve(D, Fi)
        self.e = (self.stfs.Z).dot(u+self.u_stg1)
        self.S = self.stfs.vsigma(self.e, self.prpr.e_2, self.prpr.e_0, self.prpr.e_1,
                                  self.prpr.e1, self.prpr.e0, self.prpr.e2, self.prpr.S_1, self.prpr.S1,
                                  self.prpr.Rc, self.prpr.Rt, self.prpr.E, 1)
        return u


class Convergence:
    def __init__(self, stfs, prpr, stg2):
        self.stfs = stfs
        self.prpr = prpr
        self.stg2 = stg2

    def converg(self, Fg, itr):
        if self.stg2 == 0:
            with open('fu.pickle', 'wb') as f:
                pickle.dump(itr.u, f)
        self.Ff = np.array((sum(itr.Sig*self.stfs.x*self.stfs.A), sum(itr.Sig*self.stfs.y *
                                                                      self.stfs.A), sum(itr.Sig*self.stfs.A)))
        cvr = np.hstack(
            (1000*Fg.reshape(3, 1), np.round(1000*self.Ff.reshape(3, 1), 3), itr.u.reshape(3, 1)))
        self.cvrg = pd.DataFrame(cvr, index=['Mx', 'My', 'N'], columns=[
                                 'Given', 'Found', 'u'])


class Loads:
    def __init__(self, load, prpr):
        self.load = load
        self.prpr = prpr

    def ld(self):
        if self.prpr.lim_st == '2_st' and self.prpr.dur == 'short':
            self.F = self.load.values[:, 3:6]/1000
        elif self.prpr.lim_st == '2_st' and self.prpr.dur == 'long':
            self.F = self.load.values[:, 6:9]/1000
        else:
            self.F = self.load.values[:, :3]/1000


class Result:
    def __init__(self, grade, stfs):
        self.grade = grade
        self.stfs = stfs

    def rslt1(self, itr):
        res1 = np.hstack((self.grade.reshape(self.stfs.n, 1), self.stfs.x.reshape(self.stfs.n, 1),
                          self.stfs.y.reshape(self.stfs.n, 1), itr.eps.reshape(self.stfs.n, 1), itr.Sig.reshape(self.stfs.n, 1), self.stfs.A.reshape(self.stfs.n, 1)))
        self.result1 = pd.DataFrame(res1, columns=['Grade', 'Zx', 'Zy', 'Strain', 'Stress', 'Area'])


class Interaction:
    def __init__(self, stfs, prpr):
        self.stfs = stfs
        self.prpr = prpr

    def interact(self):
        epssu = 0.025
        epsb2 = -0.0035
        Zbx = self.prpr.concr['Zx']
        Zby = self.prpr.concr['Zy']
        Zsx = self.prpr.steel['Zx']
        Zsy = self.prpr.steel['Zy']
        Zx = self.prpr.pr['Zx']
        Zy = self.prpr.pr['Zy']
        Z = np.transpose(np.array([Zx, Zy, np.ones(len(Zx))]))
        h0xmax = max(Zsx)-min(Zbx)
        rxmax = (epssu-epsb2)/h0xmax
        h0xmin = max(Zbx)-min(Zsx)
        rxmin = -(epssu-epsb2)/h0xmin
        h0ymax = max(Zsy)-min(Zby)
        rymax = (epssu-epsb2)/h0ymax
        h0ymin = max(Zby)-min(Zsy)
        rymin = -(epssu-epsb2)/h0ymin
        e0x = 0.0
        e0y = 0.0
        nd = 21
        Fi = np.zeros(0)
        r_x = np.linspace(rxmin, rxmax, nd)
        e0_x = np.absolute(min(Zbx)*r_x)
        r_y = np.linspace(rymin, rymax, nd)
        e0_y = np.absolute(min(Zby)*r_y)
        e0 = np.linspace(epsb2, epssu, nd)
        for k in (e0):
            for i in r_x:
                e0x = e0_x[np.where(r_x == i)[0][0]]
                for j in r_y:
                    e0y = e0_y[np.where(r_y == j)[0][0]]
                    u = np.array(([i, j, k+e0x+e0y]))
                    eps = Z.dot(u)
                    Sig = self.stfs.vsigma(eps, self.prpr.e_2, self.prpr.e_0, self.prpr.e_1, self.prpr.e1, self.prpr.e0,
                                           self.prpr.e2, self.prpr.S_1, self.prpr.S1, self.prpr.Rc, self.prpr.Rt, self.prpr.E, 1)
                    Ff = np.array((sum(Sig*self.stfs.x*self.stfs.A),
                                   sum(Sig*self.stfs.y * self.stfs.A), sum(Sig*self.stfs.A)))
                    Fi = np.append(Fi, Ff)
        return Fi


class Start:
    def __init__(self, fsc, prop, load, dat_file, fname):
        self.fsc = fsc
        self.prop = prop
        self.load = load
        self.dat_file = dat_file
        self.prpr = self.prop.Properties()
        self.fname = fname
        self.repname = fname+'.txt'
        self.resname = fname+'.csv'
        self.report = open(os.path.join(self.prop.rpath, self.repname), 'w')

    def strt(self, time, lst2, stg2, gb3, iF):
        start = time.time()
        print('Program =RC_frame.v1.0= designed by G. Berezin in 2020')
        print('Program =RC_frame.v1.0= designed by G. Berezin in 2020',
              file=self.report)
        print('L/C', iF, file=self.report)
        self.prpr.gb3 = gb3
        self.prpr.prp(self.dat_file)
        stfs = Stiffness(self.prpr)
        self.intrn = Interaction(stfs, self.prpr)
        self.lds = Loads(self.load, self.prpr)
        if stg2 != 0:
            with open('fu.pickle', 'rb') as f:
                u_stg1 = pickle.load(f)
        else:
            u_stg1 = np.zeros(3)
        clc = Calc(stfs, self.prpr, u_stg1)
        itr = Iteration(stfs, clc, self.prpr, self.prpr.dat)
        rsl = Result(self.prpr.grade, stfs)
        itr_b = Iteration(stfs, clc, self.prpr, self.prpr.dat)
        cnvr = Convergence(stfs, self.prpr, stg2)
        print('Results are in the files:')
        print(self.repname)
        if lst2 == 0:
            fls = self.fsc.First(self.prpr, self.lds, cnvr, itr,
                                 itr_b, clc, rsl, stfs, iF, self.dat_file, self.fname)
            fls.calc(self.report)
            print(self.resname)
        else:
            if os.path.exists(os.path.join(
                    self.prop.rpath, 'cracks.csv')):
                os.remove(os.path.join(
                    self.prop.rpath, 'cracks.csv'))
            sls = self.fsc.Second(self.prpr, self.lds, cnvr,
                                  itr, itr_b, clc, rsl, stfs, iF, self.dat_file, self.fname)
            sls.calc(self.report)
            result3 = pd.DataFrame(sls.crc.rescrack, columns=['crc', 'Abt,cm^2', 'As,cm^2', 'ds,mm', 'ls,m', 'Phi1',
                                                              'Phi2', 'Phi3', 'Psi,s', 'Sscrc,MPa', 'Ss,MPa', 'Es,MPa', 'acrc,mm'])
            result3.to_csv(os.path.join(
                self.prpr.rpath, 'cracks.csv'), sep=';', index=False)
            print('res2.csv')
            if os.path.exists(os.path.join(self.prop.rpath, 'cracks.csv')):
                print('cracks.csv')
        jt = round((time.time() - start), 4)
        print("Job time:", jt, "seconds")
        print("Job time:", jt, "seconds", file=self.report)
        self.report.close()
        return stfs, self.prpr, rsl

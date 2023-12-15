import numpy as np
import pandas as pd
import pickle
import os
pl1 = np.eye(3)
acc = 0.00000001


class Iteration:
    def __init__(self, stfs, clc, prpr, dat):
        self.clc = clc
        self.prpr = prpr
        self.stfs = stfs
        self.dat = dat

    def fcrc(self, F):
        Fcrc = pd.DataFrame(F.reshape(1, 6), columns=[
                            'Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy'])
        return Fcrc

    def d_s(self):
        n = self.dat['n'].values
        A = self.dat['A'].values/n
        ds = round(np.sqrt(4*max(A[max(self.Sig) == self.Sig])/np.pi), 3)
        return ds

    def itrn(self, Fg, iF):
        vb = np.ones((self.stfs.K, 2))
        vs = np.ones((self.stfs.ns))
        D = self.stfs.d(self.stfs.E_b, vb, self.stfs.E_s, vs, self.stfs.orientation)
        vb, Sb, Sxyb, orientation, eps1, eps2, vs, Sxys, strain, stress, u = self.clc.calc(D, Fg)
        i = 0
        du = 0.1
        self.resmes = 'Convergence is reached !'
        while du >= acc:
            i = i + 1
            if min(eps2) < -0.0035 or min(eps1) < -0.0035 or max(strain) > 0.025:
                eps1 = np.zeros(self.stfs.K)
                eps2 = np.zeros(self.stfs.K)
                strain = np.zeros(self.stfs.ns)
                self.resmes = 'Convergence L/C ' + str(iF) + ' did not reached!!!'
                print(self.resmes)
                break
            D = self.stfs.d(
                self.stfs.E_b, vb, self.stfs.E_s, vs, orientation)
            vb, Sb, Sxyb, orientation, eps1, eps2, vs, Sxys, strain, stress, u_f = self.clc.calc(
                D, Fg)
            du = np.max(abs(u - u_f))
            u = u_f
        self.eps1 = eps1
        self.eps2 = eps2
        self.Sig1 = Sb[:][:, 0]
        self.Sig2 = Sb[:][:, 1]
        self.Sxyb = Sxyb
        self.Sxys = Sxys
        self.orientation = orientation
        self.strain = strain
        self.stress = stress
        self.eps = np.append(self.eps1.reshape(
            self.stfs.K, 1), self.strain.reshape(self.stfs.ns, 1))
        self.Sig = np.append(self.Sig1.reshape(
            self.stfs.K, 1), self.stress.reshape(self.stfs.ns, 1))
        self.u = u
        self.i = i


class Stiffness:
    v = 0.0

    def __init__(self, prpr):
        self.prpr = prpr
        self.t = self.prpr.concr.A.values
        self.Zb = self.prpr.concr.Z.values
        self.K = len(self.Zb)
        self.v01 = np.ones(self.K)*self.v
        self.Eb_ = self.prpr.concr.E
        self.E_b = np.stack((self.Eb_, self.Eb_), axis=-1)
        self.orientation = np.zeros(self.K)
        self.plb = np.zeros((self.K, 3, 6))
        for i in range(0, self.K):
            pl2 = pl1 * self.Zb[i]
            self.plb[i, :, :] = np.hstack((pl1, pl2))
        self.alpha = np.radians(self.prpr.steel['alpha'].values)
        self.As = self.prpr.steel.A.values
        self.Zs = self.prpr.steel.Z.values
        self.ns = len(self.Zs)
        self.E_s = self.prpr.steel.E.values
        self.pls = np.zeros((self.ns, 3, 6))
        for i in range(0, self.ns):
            pl2 = np.eye(3) * self.Zs[i]
            self.pls[i, :, :] = np.hstack((pl1, pl2))
        self.vsigma = np.vectorize(self.prpr.sigma)

    def cQb(self, E_b):
        Qb = np.zeros((self.K, 3, 3))
        G01 = np.zeros(self.K)
        self.v10 = np.zeros(self.K)
        for i in range(0, self.K):
            if E_b[i, 0] != 0.0:
                self.v10[i] = E_b[i, 1]*self.v01[i]/E_b[i, 0]
            else:
                self.v10[i] = 0.0
            G01[i] = self.E_b[i, 0]/(2*(1+self.v10[i]))
        self.v01 = self.v10
        v = 1-self.v01*self.v10
        Qb[:, 0, 0] = E_b[:, 0]/v
        Qb[:, 0, 1] = self.v01*E_b[:, 1]/v
        Qb[:, 1, 1] = E_b[:, 1]/v
        Qb[:, 1, 0] = self.v10*E_b[:, 1]/v
        Qb[:, 2, 2] = G01
        return Qb

    def cT(self, a):
        n = len(a)
        c = np.cos(a)
        s = np.sin(a)
        T = np.zeros((n, 3, 3))
        for i in range(0, n):
            T[i, 0, :] = np.array((c[i]**2, s[i]**2, (2*c[i]*s[i])))
            T[i, 1, :] = np.array((s[i]**2, c[i]**2, (-2*c[i]*s[i])))
            T[i, 2, :] = np.array(((-c[i]*s[i]), c[i]*s[i], (c[i]**2-s[i]**2)))
        return T

    def cD(self, A, Z, T, Q):
        R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
        R_ = np.linalg.inv(R)
        T_ = np.linalg.inv(T)
        n = len(A)
        M = T_ @ Q @ R @ T @ R_
        A_ = sum(M * A.reshape(n, 1, 1))
        B_ = sum(M * Z.reshape(n, 1, 1) * A.reshape(n, 1, 1))
        D_ = sum(M * Z.reshape(n, 1, 1) ** 2 * A.reshape(n, 1, 1))
        Di = np.vstack((np.hstack((A_, B_)), np.hstack((B_, D_))))
        return Di

    def d(self, E_b, vb, E_s, vs, orientation):
        Qb = self.cQb(E_b*vb)
        T = self.cT(orientation)
        Db = self.cD(self.t, self.Zb, T, Qb)
        Qs = np.zeros((self.ns, 3, 3))
        Qs[:, 0, 0] = E_s*vs
        T = self.cT(self.alpha)
        Ds = self.cD(self.As, self.Zs, T, Qs)
        D = Db+Ds
        return D


class Calc:
    def __init__(self, stfs, prpr, u_stg1):
        self.stfs = stfs
        self.prpr = prpr
        self.u_stg1 = u_stg1

    def conc(self, u):
        v = 1-self.stfs.v01*self.stfs.v10
        K = self.stfs.K
        epsb = (self.stfs.plb @ (u+self.u_stg1)).reshape(K, 3, 1)
        exx = epsb[:, 0].reshape(K)
        eyy = epsb[:, 1].reshape(K)
        gxy = epsb[:, 2].reshape(K)
        ee1 = exx + eyy
        ee2 = exx - eyy
        emax = (ee1/2 + np.sqrt((ee2 / 2) ** 2 + (gxy/2) ** 2))
        emin = (ee1/2 - np.sqrt((ee2 / 2) ** 2 + (gxy/2) ** 2))
        eps1 = (emax + self.stfs.v01*emin)/v
        eps2 = (self.stfs.v10*emax + emin)/v
        orientation = 0.5 * np.arctan2(gxy, ee2)
        kRb = np.ones(K)
        for i in range(0, K):
            if eps1[i] > 0.002 and eps2[i] < 0:
                kRb[i] = 1.0 / (0.8 + 100 * eps1[i])
            else:
                kRb[i] = 1.0
        Sb = np.vstack((self.stfs.vsigma(eps1, self.prpr.ce_2, self.prpr.ce_0, self.prpr.ce_1, self.prpr.ce1, self.prpr.ce0, self.prpr.ce2, self.prpr.cS_1, self.prpr.cS1, self.prpr.cRc*kRb, self.prpr.cRt, self.prpr.cE, 1),
                        self.stfs.vsigma(eps2, self.prpr.ce_2, self.prpr.ce_0, self.prpr.ce_1, self.prpr.ce1, self.prpr.ce0, self.prpr.ce2, self.prpr.cS_1, self.prpr.cS1, self.prpr.cRc*kRb, self.prpr.cRt, self.prpr.cE, kRb[i]))).transpose().reshape(K, 2, 1)
        vb = self.v_b(K, Sb, eps1, eps2)
        Sxyb = self.sxyb(K, orientation, Sb)
        return vb, Sb, Sxyb, orientation, eps1, eps2

    def sxyb(self, K, orientation, Sb):
        c = np.cos(orientation)
        s = np.sin(orientation)
        Sxyb = np.zeros((K, 3, 1))
        for i in range(0, K):
            Cb = np.array(
                [[c[i] ** 2, s[i] ** 2], [s[i] ** 2, c[i] ** 2], [s[i] * c[i], -s[i] * c[i]]])
            Sxyb[i][:, :] = Cb @ Sb[i][:, :]
        return Sxyb

    def v_b(self, K, Sb, eps1, eps2):
        vb = np.ones((K, 2))
        for i in range(0, K):
            if eps1[i] != 0:
                vb[i, 0] = Sb[i][0] / self.stfs.E_b[i][0]/eps1[i]
            else:
                vb[i, 0] = 1.0
            if eps2[i] != 0:
                vb[i, 1] = Sb[i][1] / self.stfs.E_b[i][1] / eps2[i]
            else:
                vb[i, 1] = 1.0
        return vb

    def reb(self, u):
        ns = self.stfs.ns
        epss = (self.stfs.pls @ (u+self.u_stg1)).reshape(ns, 3, 1)
        c = np.cos(self.stfs.alpha)
        s = np.sin(self.stfs.alpha)
        strain = np.zeros((ns))
        for i in range(0, ns):
            dc = np.array([c[i] ** 2, s[i] ** 2, 2 * s[i] * c[i]])
            strain[i] = dc @ epss[i]
        stress = self.stfs.vsigma(
            strain, self.prpr.re_2, self.prpr.re_0, self.prpr.re_1, self.prpr.re1, self.prpr.re0, self.prpr.re2, self.prpr.rS_1, self.prpr.rS1, self.prpr.rRc, self.prpr.rRt, self.prpr.rE, 1)
        vs = self.v_s(ns, strain, stress)
        Sxys = self.sxys(ns, c, s, stress)
        return vs, Sxys, strain, stress

    def sxys(self, ns, c, s, stress):
        Sxys = np.zeros((ns, 3, 1))
        for i in range(0, ns):
            Cs = np.array([[c[i] ** 2, s[i] ** 2], [s[i] ** 2,
                                                    c[i] ** 2], [s[i] * c[i], -s[i] * c[i]]])
            STs = np.hstack([stress[i], 0])
            Sxys[i][:, :] = (Cs @ STs).reshape(3, 1)
        return Sxys

    def v_s(self, ns, strain, stress):
        vs = np.ones((ns))
        for i in range(0, ns):
            if strain[i] != 0:
                vs[i] = stress[i] / self.stfs.E_s[i] / strain[i]
            else:
                vs[i] = 1.0
        return vs

    def calc(self, D, F):
        try:
            u = np.linalg.solve(D, F)
            self.mes = 'ok'
        except np.linalg.LinAlgError as var1:
            self.mes = var1
            D = np.eye(6)
            Fi = np.array(([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
            u = np.linalg.solve(D, Fi)
        vb, Sb, Sxyb, orientation, eps1, eps2 = self.conc(u)
        vs, Sxys, strain, stress = self.reb(u)
        return vb, Sb, Sxyb, orientation, eps1, eps2, vs, Sxys, strain, stress, u


class Convergence:
    def __init__(self, stfs, prpr, stg2):
        self.stfs = stfs
        self.prpr = prpr
        self.stg2 = stg2

    def converg(self, Fg, itr):
        Z = self.prpr.pr['Z'].values
        Sxy = np.vstack((itr.Sxyb, itr.Sxys))
        A = self.prpr.pr['A'].values
        if self.stg2 == 0:
            with open('su.pickle', 'wb') as f:
                pickle.dump(itr.u, f)
        F = np.zeros(6)
        n = len(Z)
        for i in range(0, n):
            F[0] = F[0] + Sxy[i][0] * A[i]
            F[1] = F[1] + Sxy[i][1] * A[i]
            F[2] = F[2] + Sxy[i][2] * A[i]
            F[3] = F[3] + Sxy[i][0] * A[i] * Z[i]
            F[4] = F[4] + Sxy[i][1] * A[i] * Z[i]
            F[5] = F[5] + Sxy[i][2] * A[i] * Z[i]
        cvr = np.hstack((Fg.reshape(6, 1), np.round(
            F.reshape(6, 1), 4), itr.u.reshape(6, 1)))
        self.cvrg = pd.DataFrame(cvr, index=[
                                 'Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy'], columns=['Given', 'Found', 'u'])


class Loads:
    def __init__(self, load, prpr):
        self.load = load
        self.prpr = prpr

    def ld(self):
        self.lim_st = self.prpr.lim_st
        self.dur = self.prpr.dur
        if self.lim_st == '2_st' and self.dur == 'short':
            self.F = self.load.values[:, 6:12]/1000
        elif self.lim_st == '2_st' and self.dur == 'long':
            self.F = self.load.values[:, 12:18]/1000
        else:
            self.F = self.load.values[:, :6]/1000


class Result:
    def __init__(self, stfs, prpr):
        self.stfs = stfs
        self.prpr = prpr

    def rslt1(self, itr):
        ang_c = np.degrees(itr.orientation)
        ang_s = np.degrees(self.stfs.alpha)
        res1_b = np.hstack((self.prpr.grade_b.reshape(self.stfs.K, 1),
                            self.stfs.Zb.reshape(self.stfs.K, 1),
                            itr.eps1.reshape(self.stfs.K, 1),
                            itr.eps2.reshape(self.stfs.K, 1),
                            itr.Sig1.reshape(self.stfs.K, 1),
                            itr.Sig2.reshape(self.stfs.K, 1),
                            np.round(ang_c.reshape(self.stfs.K, 1), 3),
                            self.stfs.t.reshape(self.stfs.K, 1)))
        res1_s = np.hstack((self.prpr.grade_s.reshape(self.stfs.ns, 1),
                            self.stfs.Zs.reshape(self.stfs.ns, 1),
                            itr.strain.reshape(self.stfs.ns, 1),
                            np.zeros(self.stfs.ns).reshape(self.stfs.ns, 1),
                            itr.stress.reshape(self.stfs.ns, 1),
                            np.zeros(self.stfs.ns).reshape(self.stfs.ns, 1),
                            np.round(ang_s.reshape(self.stfs.ns, 1), 3),
                            self.stfs.As.reshape(self.stfs.ns, 1)))
        res1 = np.vstack((res1_b, res1_s))

        self.result1 = pd.DataFrame(res1, columns=[
                                    'Grade', 'Z', 'Strain1', 'Strain2', 'Stress1', 'Stress2', 'Angle', 'Area'])


class Start:
    def __init__(self, fsc, prop, load, dat_file, v, fname):
        self.fsc = fsc
        self.prop = prop
        self.load = load
        self.dat_file = dat_file
        self.v = v
        self.fname = fname
        self.repname = fname+'.txt'
        self.resname = fname+'.csv'

    def strt(self, time, lst2, stg2, gb3, iF):
        start = time.time()
        self.prpr = self.prop.Properties()
        self.report = open(os.path.join(self.prop.dpath, self.repname), 'w')
        Stiffness.v = self.v
        print('Program =RC_slab.v1.0= designed by G. Berezin in 2020')
        print('Program =RC_slab.v1.0= designed by G. Berezin in 2020',
              file=self.report)
        print('L/C', iF, file=self.report)
        self.lds = Loads(self.load, self.prpr)
        self.prpr.gb3 = gb3
        self.prpr.prp(self.dat_file)
        stfs = Stiffness(self.prpr)
        if stg2 != 0:
            with open('su.pickle', 'rb') as f:
                u_stg1 = pickle.load(f)
        else:
            u_stg1 = np.zeros(6)
        clc = Calc(stfs, self.prpr, u_stg1)
        itr = Iteration(stfs, clc, self.prpr, self.prpr.dat)
        rsl = Result(stfs, self.prpr)
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
                    self.prop.dpath, 'cracks.csv')):
                os.remove(os.path.join(
                    self.prop.dpath, 'cracks.csv'))
            sls = self.fsc.Second(self.prpr, self.lds, cnvr,
                                  itr, itr_b, clc, rsl, stfs, iF, self.dat_file)
            sls.calc(self.report)
            result3 = pd.DataFrame(sls.crc.rescrack, columns=['crc', 'Abt,cm^2', 'As,cm^2', 'ds,mm', 'ls,m', 'Phi1',
                                                              'Phi2', 'Phi3', 'Psi,s', 'Sscrc,MPa', 'Ss,MPa', 'Es,MPa', 'acrc,mm'])
            result3.to_csv(os.path.join(
                self.prpr.dpath, 'cracks.csv'), sep=';', index=False)
            print('res2.csv')
            if os.path.exists(os.path.join(self.prop.dpath, 'cracks.csv')):
                print('cracks.csv')
        jt = round((time.time() - start), 4)
        print("Job time:", jt, "seconds")
        print("Job time:", jt, "seconds", file=self.report)
        self.report.close()
        return stfs, self.prpr, rsl

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import spatial as sp_spatial
import mpl_toolkits.mplot3d as a3

class Plots:
    @staticmethod
    def strainstress(stfs, prpr, mat):
        i = np.where(mat == prpr.pr['T'])[0][0]
        x = [prpr.e_2[i], prpr.e_0[i], prpr.e_1[i],
             0.0, prpr.e1[i], prpr.e0[i], prpr.e2[i]]
        y = stfs.vsigma(x, prpr.e_2[i], prpr.e_0[i], prpr.e_1[i], prpr.e1[i], prpr.e0[i],
                        prpr.e2[i], prpr.S_1[i], prpr.S1[i], prpr.Rc[i], prpr.Rt[i], prpr.E[i], 1)
        fig, ax = plt.subplots(num=mat, dpi=150)
        ax.plot(x, y)
        ax.set_xlabel('Strain')
        ax.set_ylabel('Stress, MPa')
        plt.title('Strain-Stress')
        ax.annotate(round(y[1], 3), (x[1], y[1]), ha='center')
        ax.annotate(round(y[2], 3), (x[2], y[2]), ha='right')
        ax.annotate(round(y[3], 3), (x[3], y[3]), ha='right')
        ax.annotate(round(y[4], 3), (x[4], y[4]), ha='right')
        ax.annotate(round(y[5], 3), (x[5], y[5]), ha='left')
        ax.scatter(x, y, c='red', alpha=0.5)
        plt.show()

    @staticmethod
    def f_strain(rsl, stfs):
        xc = rsl.result1.Zx[stfs.prpr.pr['T'] == 'concrete']
        yc = rsl.result1.Zy[stfs.prpr.pr['T'] == 'concrete']
        zc = rsl.result1.Strain[stfs.prpr.pr['T'] == 'concrete']
        epsmin = min(zc)
        x_min = np.array(xc[zc == epsmin][:1])[0]
        y_min = np.array(yc[zc == epsmin][:1])[0]
        x1 = np.linspace(xc.min(), xc.max(), len(xc))
        y1 = np.linspace(yc.min(), yc.max(), len(yc))
        x2, y2 = np.meshgrid(x1, y1)
        z2 = griddata((xc, yc), zc, (x1[None, :], y1[:, None]), method='linear')
        fig, ax = plt.subplots()
        ax.set_xlabel('Zbx')
        ax.set_ylabel('Zby')
        ax.set_aspect('equal')
        plt.subplots_adjust(left=0.185, right=0.815, bottom=0.1, top=0.9)
        plt.title('Strain')
        cont = ax.contourf(x2, y2, z2, 20, alpha=0.6, cmap="rainbow")
        xr = rsl.result1.Zx[stfs.prpr.pr['T'] == 'rebar'].values
        yr = rsl.result1.Zy[stfs.prpr.pr['T'] == 'rebar'].values
        zr = rsl.result1.Strain[stfs.prpr.pr['T'] == 'rebar'].values
        As = stfs.A[stfs.prpr.pr['T'] == 'rebar']
        ds = np.sqrt(4*As/np.pi)
        ax.scatter(xr, yr, s=ds * 8000, c='green', alpha=0.5)
        ax.scatter(x_min, y_min, s=500, c='white')
        ax.annotate(round(epsmin, 5), (x_min, y_min), size=10, xytext=(
            0, 10), ha='right', color='red', textcoords='offset points')
        for i, txt in enumerate(zr):
            ax.annotate(round(txt, 5), (xr[i], yr[i]), size=10, xytext=(
                0, 0), ha='left', textcoords='offset points')
        plt.show()
        
    @staticmethod
    def f_stress(rsl, stfs):
        xc = rsl.result1.Zx[stfs.prpr.pr['T'] == 'concrete'].values
        yc = rsl.result1.Zy[stfs.prpr.pr['T'] == 'concrete'].values
        zc = rsl.result1.Stress[stfs.prpr.pr['T'] == 'concrete'].values
        stressmin = min(zc)
        x_min = np.array(xc[zc == stressmin][:1])[0]
        y_min = np.array(yc[zc == stressmin][:1])[0]
        x1 = np.linspace(xc.min(), xc.max(), len(xc))
        y1 = np.linspace(yc.min(), yc.max(), len(yc))
        x2, y2 = np.meshgrid(x1, y1)
        z2 = griddata((xc, yc), zc, (x1[None, :], y1[:, None]), method='linear')
        fig, ax = plt.subplots()
        ax.set_xlabel('Zbx')
        ax.set_ylabel('Zby')
        ax.set_aspect('equal')
        plt.subplots_adjust(left=0.185, right=0.815, bottom=0.1, top=0.9)
        plt.title('Stress, MPa')
        cont = ax.contourf(x2, y2, z2, 20, alpha=0.6, cmap="rainbow")
        xr = rsl.result1.Zx[stfs.prpr.pr['T'] == 'rebar'].values
        yr = rsl.result1.Zy[stfs.prpr.pr['T'] == 'rebar'].values
        zr = rsl.result1.Stress[rsl.stfs.prpr.pr['T'] == 'rebar'].values
        As = stfs.A[stfs.prpr.pr['T'] == 'rebar']
        ds = np.sqrt(4*As/np.pi)
        ax.scatter(xr, yr, s=ds * 8000, c='green', alpha=0.5)
        ax.scatter(x_min, y_min, s=500, c='white')
        ax.annotate(round(stressmin, 2), (x_min, y_min), size=10, xytext=(
            0, 10), ha='right', color='red', textcoords='offset points')
        for i, txt in enumerate(zr):
            ax.annotate(round(txt, 2), (xr[i], yr[i]), size=10, xytext=(
                0, 0), ha='left', textcoords='offset points')
        plt.show()

    @staticmethod
    def f_stress3D(rsl, stfs):
        xc = rsl.result1.Zx[stfs.prpr.pr['T'] == 'concrete'].values
        yc = rsl.result1.Zy[stfs.prpr.pr['T'] == 'concrete'].values
        zc = rsl.result1.Stress[stfs.prpr.pr['T'] == 'concrete'].values
        x1 = np.linspace(xc.min(), xc.max(), len(xc))
        y1 = np.linspace(yc.min(), yc.max(), len(yc))
        x2, y2 = np.meshgrid(x1, y1)
        fig = plt.figure(num='Stress3D', dpi=150)
        ax = plt.axes(projection='3d')
        z2 = griddata((xc, yc), zc, (x2, y2), method='cubic')
        surf = ax.plot_surface(x2, y2, z2, cmap=cm.jet, alpha=0.8, linewidth=1)
        plt.title('Concrete Stress, MPa')
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()

    @staticmethod
    def chart_I(Fi, Force):
        Mx = Fi[0::3]
        My = Fi[1::3]
        N = Fi[2::3]
        n = len(Mx)
        points = np.hstack(
            (Mx.reshape(n, 1), My.reshape(n, 1), N.reshape(n, 1)))
        fig = plt.figure(num='Interaction Surface', dpi=150)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(0, 0, N.max(), c='green', alpha=0.5)
        ax.text(0.0, 0.0, N.max(),  '%s' % (str(round(N.max(), 3))),
                size=10, zorder=1, color='k', ha='left', va='bottom')
        ax.scatter(0, 0, N.min(), c='green', alpha=0.5)
        ax.text(0.0, 0.0, N.min(),  '%s' % (str(round(N.min(), 3))),
                size=10, zorder=1, color='k', ha='left', va='top')
        ax.scatter(Mx.min(), 0.0, N[np.where(
            Mx == Mx.min())[0][0]], c='green', alpha=0.5)
        ax.text(Mx.min(), 0.0, N[np.where(Mx == Mx.min())[0][0]],  '%s' % (
            str(round(Mx.min(), 3))), size=10, zorder=1, color='k', ha='right',
            va='center')
        ax.scatter(Mx.max(), 0.0, N[np.where(
            Mx == Mx.max())[0][0]], c='green', alpha=0.5)
        ax.text(Mx.max(), 0.0, N[np.where(Mx == Mx.max())[0][0]],  '%s' % (
            str(round(Mx.max(), 3))), size=10, zorder=1, color='k', ha='left',
            va='center')
        ax.scatter(0.0, My.min(), N[np.where(
            My == My.min())[0][0]], c='green', alpha=0.5)
        ax.text(0.0, My.min(), N[np.where(My == My.min())[0][0]],  '%s' % (
            str(round(My.min(), 3))), size=10, zorder=1, color='k', ha='right',
            va='center')
        ax.scatter(0.0, My.max(), N[np.where(
            My == My.max())[0][0]], c='green', alpha=0.5)
        ax.text(0.0, My.max(), N[np.where(My == My.max())[0][0]],  '%s' % (
            str(round(My.max(), 3))), size=10, zorder=1, color='k', ha='left',
            va='center')
        M_x = Force[0]
        M_y = Force[1]
        N_ = Force[2]
        ax.scatter(M_x, M_y, N_, c='red', alpha=1.0)
        ax.text(M_x, M_y, N_,  '%s' % (str(Force)), size=12, zorder=1,
                color='red', ha='center', va='bottom')
        ax.set_xlabel('Mx, MN*m')
        ax.set_ylabel('My, MN*m')
        ax.set_zlabel('N, MN')
        fig.tight_layout()
        hull = sp_spatial.ConvexHull(points)
        indices = hull.simplices
        faces = points[indices]
        for f in faces:
            face = a3.art3d.Poly3DCollection([f])
            face.set_edgecolor('k')
            face.set_linewidths(0.1)
            face.set_alpha(0.5)
            ax.add_collection3d(face)
        plt.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.9)
        plt.show()

    @staticmethod
    def sc_strain(rsl, stfs):
        ax = plt.gca()
        df = rsl.result1[stfs.prpr.pr['T'] == 'concrete']
        df.plot(kind='line', x='Z', y='Strain1', color='green', ax=ax)
        df.plot(kind='line', x='Z', y='Strain2', color='red', ax=ax)
        plt.title('Strain')
        ax.set_xlabel('Z, m')
        ax.set_ylabel('Strain')
        plt.show()

    @staticmethod
    def sc_stress(rsl, stfs):
        fig = plt.figure()
        ax = plt.gca()
        df = rsl.result1[stfs.prpr.pr['T'] == 'concrete']
        df.plot(kind='line', x='Z', y='Stress1', color='green', ax=ax)
        df.plot(kind='line', x='Z', y='Stress2', color='red', ax=ax)
        plt.title('Stress, MPa')
        ax.set_xlabel('Z, m')
        ax.set_ylabel('Stress, MPa')
        plt.show()

    @staticmethod
    def strains(result):
        ax = plt.gca()
        df = result
        df.plot(kind='line', x=None, y='eps_min', color='green', ax=ax)
        df.plot(kind='line', x=None, y='eps_max', color='red', ax=ax)
        plt.title('Strain')
        ax.set_xlabel('L/C')
        ax.set_ylabel('Strain')
        plt.show()

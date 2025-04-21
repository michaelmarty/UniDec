import numpy as np
import matplotlib.pyplot as plt


class Ion:
    def __init__(self, m=150, z=1, dt=1e-9):
        self.dt = dt
        self.m = m
        self.da2kg = 1.66054e-27
        self.mass = self.m * self.da2kg
        self.z = z
        self.mz = self.m / self.z
        self.e = 1.602176634e-19
        self.ke = 8.99e9
        self.q = self.z * self.e
        self.pos = np.array([0., 0., 0.])
        self.vel = np.array([0., 0., 0.])
        self.acc = np.array([0., 0., 0.])
        self.traj = np.array([self.pos])
        self.times = np.array([0])
        self.alive = True

    def applyA(self):
        if self.alive:
            dv = self.acc * self.dt
            self.vel = self.vel + dv
            self.pos = self.pos + self.vel * self.dt
            self.traj = np.append(self.traj, [self.pos], axis=0)
            t = self.times[-1] + self.dt
            self.times = np.append(self.times, t)
        return self.pos

    def applyE(self, Efield):
        self.acc = self.q * Efield / self.mass
        return self.applyA()

    def applyPS(self, Ions):
        const = self.ke * self.q / self.mass
        p1 = self.pos
        a = np.array([0.0, 0.0, 0.0])
        for i in Ions:
            if i.alive:
                p2 = i.pos
                q2 = i.q
                dist = np.linalg.norm(p1 - p2)
                if dist != 0:
                    aval = const * q2 / (dist ** 2)
                    a += aval
        self.acc = a
        return self.applyA()


class Ions:
    def __init__(self, mlist, zlist, dt=1e-9):
        self.ions = [Ion(m, zlist[i], dt) for i, m in enumerate(mlist)]

    def MonteCarlo(self, dxyz, dv):
        for ion in self.ions:
            random_pos = np.random.normal(scale=dxyz, size=3)
            ion.pos += random_pos
            ion.traj = np.array([ion.pos])
            random_v = np.random.normal(scale=dv, size=3)
            ion.vel += random_v

    def ToF(self):
        Efield1 = np.array([1e6, 100000, 0])
        Efield2 = np.array([-1e6, 0, 0])
        drift = np.array([0, 0, 0])
        for i in np.arange(30000):
            for ion in self.ions:
                if ion.alive:
                    ion.applyPS(self.ions)
                    if ion.pos[0] < 0.01:
                        Efield = Efield1
                    elif ion.pos[0] > 1:
                        Efield = Efield2
                    else:
                        Efield = drift
                    if ion.pos[1] > 0.01 and ion.pos[0] < 0.01:
                        ion.alive = False
                        break
                    ion.applyE(Efield)


class Efield:
    def __init__(self):
        e0 = 8.85e-12  # C^2/Nm^2


if __name__ == "__main__":
    mlist = [150, 150, 150, 150]
    zlist = [1, 1, 1, 1]

    ions = Ions(mlist, zlist)
    ions.MonteCarlo(0.001, 1)
    ions.ToF()

    for ion in ions.ions:
        plt.subplot(131)
        plt.plot(ion.times, ion.traj[:, 0], marker="o")
        plt.subplot(132)
        plt.plot(ion.times, ion.traj[:, 1], marker="o")
        plt.subplot(133)
        plt.plot(ion.traj[:, 0], ion.traj[:, 1], marker="o")
    plt.show()

    exit()

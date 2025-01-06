import matplotlib.pyplot as plt
import numpy as np

def acceleration(rVec, M, G=6.674e-11):
    r = np.linalg.norm(rVec)

    if r > 0:
        return - G * M * rVec / r**3
    else:
        print('Error: Division by 0, return 0')
        return np.zeros_like(rVec)

class Celestial(object):
    def __init__(self, startPosition, startVelocity, mass, N):
        self.position = np.zeros((N, 2), dtype=object)
        self.position[0] = np.array(startPosition, dtype=float)
        self.velocity = np.zeros((N, 2), dtype=object)
        self.velocity[0] = np.array(startVelocity, dtype=float)
        self.mass = mass

    def update(self, celestial_bodies, i, dt):
        self.velocity[i] = self.velocity[i-1]
        for other in celestial_bodies:
            self.velocity[i] = self.velocity[i] + dt * acceleration(rVec=(self.position[i-1] - other.position[i-1]), M=other.mass)

        self.position[i] = self.position[i-1] + dt * self.velocity[i]

class SolarSystem(object):
    def __init__(self, N, dt):
        self.Earth = Celestial(startPosition=[149.6e9, 0], startVelocity=[0, 29.8e3], mass=0.330e24, N=N)
        self.Sun = Celestial(startPosition=[0, 0], startVelocity=[0, 0], mass=1.989e30, N=N)
        self.Jupiter = Celestial(startPosition=[778.5e9, 0], startVelocity=[0, 13.1e3], mass=1898e24, N=N)

        self.N = N
        self.dt = dt

    def simulate(self):
        for i in range(1, self.N):
            self.Earth.update([self.Sun, self.Jupiter], i, self.dt)
            self.Jupiter.update([self.Sun, self.Earth], i, self.dt)

    def savePlot(self):
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(0, 0, marker='o', color='tab:orange', label='Sun')
        ax.plot(self.Earth.position[:, 0], self.Earth.position[:, 1], linestyle='-', color='tab:blue', label='Earth')
        ax.plot(self.Jupiter.position[:, 0], self.Jupiter.position[:, 1], linestyle='-', color='tab:red', label='Jupiter')
        ax.plot(self.Earth.position[0, 0], self.Earth.position[0, 1], linestyle=None, marker='o', color='tab:blue')
        ax.plot(self.Jupiter.position[0, 0], self.Jupiter.position[0, 1], linestyle=None, marker='o', color='tab:red')
        ax.plot(self.Earth.position[-1, 0], self.Earth.position[-1, 1], linestyle=None, marker='h', color='tab:blue')
        ax.plot(self.Jupiter.position[-1, 0], self.Jupiter.position[-1, 1], linestyle=None, marker='h', color='tab:red')

        ax.set(xlim=(-1e12, 1e12), 
               ylim=(-1e12, 1e12), 
               xlabel='x (m)',
               ylabel='y (m)')
        ax.set_title('Simulation length: 11.5 years')
        ax.legend(loc=2)

        fig.savefig('Figures/simple_solar_system.png', dpi=600)

def main():
    sim = SolarSystem(N=int(365*11.5), dt=60 * 60 * 24)
    sim.simulate()
    sim.savePlot()

if __name__ == "__main__":
    main()
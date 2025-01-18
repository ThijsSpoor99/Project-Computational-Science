import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
from joblib import Parallel, delayed

def acceleration(rVec, M, G=6.674e-11):
    r = np.linalg.norm(rVec)
    if r > 0:
        return - G * M * rVec / r**3
    else:
        print('Error: Division by 0, return 0')
        return np.zeros_like(rVec)

class Celestial(object):
    def __init__(self, startPosition, startVelocity, mass, N):
        self.position = np.zeros((N, 3), dtype=float)
        self.position[0] = np.array(startPosition, dtype=float)
        self.velocity = np.zeros((N, 3), dtype=float)
        self.velocity[0] = np.array(startVelocity, dtype=float)
        self.mass = mass

    def update(self, celestial_bodies, i, dt):
        self.velocity[i] = self.velocity[i-1]
        for other in celestial_bodies:
            self.velocity[i] = self.velocity[i] + dt * acceleration(
                rVec=(self.position[i-1] - other.position[i-1]), M=other.mass
            )
        self.position[i] = self.position[i-1] + dt * self.velocity[i]

class SolarSystem:
    def __init__(self, N, dt):
        self.planets = {}
        startPositionsPlanets = pd.read_csv("Data/planetData.csv", header=0, index_col=0)
        for planet, pos in startPositionsPlanets.iterrows():
            self.planets[planet] = Celestial(
                startPosition=pos.loc[['x', 'y', 'z']].values,
                startVelocity=pos.loc[['Vx', 'Vy', 'Vz']].values,
                mass=pos.loc['M'],
                N=N
            )

            self.planetColors = {"Sun": "yellow",
                                 "Mercury": "gray",
                                 "Venus": "khaki",
                                 "Earth": "blue",
                                 "Mars": "red",
                                 "Jupiter": "darkorange",
                                 "Saturn": "gold",
                                 "Uranus": "lightskyblue",
                                 "Neptune": "royalblue"}

        self.asteroids = {}
        startPositionsAsteroids = pd.read_csv("Data/CentaursCartesian.csv", header=0, index_col=0).head(50)
        for asteroid, pos in startPositionsAsteroids.iterrows():
            self.asteroids[asteroid] = Celestial(
                startPosition=pos.loc[['x', 'y', 'z']].values,
                startVelocity=pos.loc[['Vx', 'Vy', 'Vz']].values,
                mass=1,  # Massless assumption
                N=N
            )

        self.planets["Sun"] = Celestial(
            startPosition=[0, 0, 0],
            startVelocity=[0, 0, 0],
            mass=1.989e30,
            N=N
        )
        self.N = N
        self.dt = dt
    
    def simulate(self):
        for i in range(1, self.N):
            # update planets
            for name, body in self.planets.items():
                other_planets = [obj for planet, obj in self.planets.items() if planet != name]
                body.update(other_planets, i, self.dt)
            
            # update asteroids
            for name, body in self.asteroids.items():
                body.update(self.planets.values(), i, self.dt)

    def animate(self):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1, 1, 1))

        # creating axes for planets and asteroids
        plots = {}
        trails = {}
        for name in self.planets.keys():
            plots[name], = ax.plot([], [], [], 'o', label=name, markersize=5, c=self.planetColors[name])
            trails[name], = ax.plot([], [], [], linestyle='-', c=self.planetColors[name])

        for name in self.asteroids.keys():
            plots[name], = ax.plot([], [], [], 'o', markersize=2, c="dimgray")
            trails[name], = ax.plot([], [], [], linestyle=':', c="darkgray")

        def init():
            """parameters for 3d plot"""
            ax.set_xlim([-5e12, 5e12])
            ax.set_ylim([-5e12, 5e12])
            ax.set_zlim([-5e12, 5e12])
            ax.set_xlabel('x (m)')
            ax.set_ylabel('y (m)')
            ax.set_zlabel('z (m)')
            ax.set_title('Live 3D Solar System Simulation')
            ax.legend()
            return [*plots.values(), *trails.values()]

        def update(frame):
            #update planets and their trails
            for name, body in self.planets.items():
                plots[name].set_data([body.position[frame, 0]], [body.position[frame, 1]])
                plots[name].set_3d_properties([body.position[frame, 2]])
                trails[name].set_data(body.position[:frame, 0], body.position[:frame, 1])
                trails[name].set_3d_properties(body.position[:frame, 2])
            
            #update asteroids and their (shortened) trails
            for name, body in self.asteroids.items():
                plots[name].set_data([body.position[frame, 0]], [body.position[frame, 1]])
                plots[name].set_3d_properties([body.position[frame, 2]])

                start_idx = max(0, frame - 3000) #3000 is maximum trail length
                trails[name].set_data(body.position[start_idx:frame, 0], body.position[start_idx:frame, 1])
                trails[name].set_3d_properties(body.position[start_idx:frame, 2])
            return [*plots.values(), *trails.values()]

        ani = FuncAnimation(fig, update, frames=range(0, self.N, 200), init_func=init, blit=False, interval=100)
        ani.save("Figures/solar_system_3D.gif", writer=PillowWriter(fps=30))

def main():
    sim = SolarSystem(N=int(365 * 100), dt=60 * 60 * 24)  # simulate for 100 years
    sim.simulate()
    sim.animate()

if __name__ == "__main__":
    main()
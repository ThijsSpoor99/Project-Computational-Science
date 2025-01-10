import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter

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

class SolarSystem(object):
    def __init__(self, N, dt):
        startPositions = pd.read_csv("Data/planetData.csv", header=0, index_col=0)
        for planet, pos in startPositions.iterrows():
            setattr(self, planet,  Celestial(
                startPosition= pos.loc[['x','y','z']].values,
                startVelocity= pos.loc[['Vx','Vy','Vz']].values,
                mass=pos.loc[['M']].values,
                N=N
            ))

        self.Sun = Celestial(
            startPosition=[0, 0, 0],
            startVelocity=[0, 0, 0],
            mass=1.989e30,
            N=N
        )

        self.N = N
        self.dt = dt

    def simulate(self):
        for i in range(1, self.N):
            self.Earth.update([self.Sun, self.Jupiter, self.Mars, self.Saturn, self.Neptune], i, self.dt)
            self.Jupiter.update([self.Sun, self.Earth, self.Mars, self.Saturn, self.Uranus, self.Neptune], i, self.dt)
            self.Mars.update([self.Sun, self.Earth, self.Jupiter, self.Saturn, self.Uranus, self.Neptune], i, self.dt)
            self.Saturn.update([self.Sun, self.Earth, self.Jupiter, self.Mars, self.Uranus, self.Neptune], i, self.dt)
            self.Uranus.update([self.Sun, self.Earth, self.Jupiter, self.Mars, self.Saturn, self.Neptune], i, self.dt)
            self.Neptune.update([self.Sun, self.Earth, self.Jupiter, self.Mars, self.Saturn, self.Uranus], i, self.dt)

    def animate(self):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1, 1, 1))

        sun_plot, = ax.plot([], [], [], 'o', color='tab:orange', label='Sun', markersize=10)
        earth_plot, = ax.plot([], [], [], 'o', color='tab:blue', label='Earth', markersize=5)
        jupiter_plot, = ax.plot([], [], [], 'o', color='tab:red', label='Jupiter', markersize=8)
        mars_plot, = ax.plot([], [], [], 'o', color='tab:brown', label='Mars', markersize=5)
        saturn_plot, = ax.plot([], [], [], 'o', color='tab:olive', label='Saturn', markersize=7)
        uranus_plot, = ax.plot([], [], [], 'o', color='tab:cyan', label='Uranus', markersize=6)
        neptune_plot, = ax.plot([], [], [], 'o', color='tab:purple', label='Neptune', markersize=6)
        
        earth_trail, = ax.plot([], [], [], linestyle='-', color='tab:blue')
        jupiter_trail, = ax.plot([], [], [], linestyle='-', color='tab:red')
        mars_trail, = ax.plot([], [], [], linestyle='-', color='tab:brown')
        mars_trail, = ax.plot([], [], [], linestyle='-', color='tab:brown')
        saturn_trail, = ax.plot([], [], [], linestyle='-', color='tab:olive')
        uranus_trail, = ax.plot([], [], [], linestyle='-', color='tab:cyan')
        neptune_trail, = ax.plot([], [], [], linestyle='-', color='tab:purple')

        def init():
            ax.set_xlim([-5e12, 5e12])
            ax.set_ylim([-5e12, 5e12])
            ax.set_zlim([-0.8e12, 0.8e12])
            ax.set_xlabel('x (m)')
            ax.set_ylabel('y (m)')
            ax.set_zlabel('z (m)')
            ax.set_title('Live 3D Solar System Simulation')
            ax.legend()
            return sun_plot, earth_plot, jupiter_plot, mars_plot, saturn_plot, uranus_plot, neptune_plot, earth_trail, jupiter_trail, mars_trail, saturn_trail, uranus_trail, neptune_trail

        def update(frame):
            # update positions
            sun_plot.set_data([0], [0])
            sun_plot.set_3d_properties([0])

            earth_plot.set_data([self.Earth.position[frame, 0]], [self.Earth.position[frame, 1]])
            earth_plot.set_3d_properties([self.Earth.position[frame, 2]])

            jupiter_plot.set_data([self.Jupiter.position[frame, 0]], [self.Jupiter.position[frame, 1]])
            jupiter_plot.set_3d_properties([self.Jupiter.position[frame, 2]])

            mars_plot.set_data([self.Mars.position[frame, 0]], [self.Mars.position[frame, 1]])
            mars_plot.set_3d_properties([self.Mars.position[frame, 2]])

            saturn_plot.set_data([self.Saturn.position[frame, 0]], [self.Saturn.position[frame, 1]])
            saturn_plot.set_3d_properties([self.Saturn.position[frame, 2]])

            uranus_plot.set_data([self.Uranus.position[frame, 0]], [self.Uranus.position[frame, 1]])
            uranus_plot.set_3d_properties([self.Uranus.position[frame, 2]])

            neptune_plot.set_data([self.Neptune.position[frame, 0]], [self.Neptune.position[frame, 1]])
            neptune_plot.set_3d_properties([self.Neptune.position[frame, 2]])

            # update trails
            earth_trail.set_data(self.Earth.position[:frame, 0], self.Earth.position[:frame, 1])
            earth_trail.set_3d_properties(self.Earth.position[:frame, 2])

            jupiter_trail.set_data(self.Jupiter.position[:frame, 0], self.Jupiter.position[:frame, 1])
            jupiter_trail.set_3d_properties(self.Jupiter.position[:frame, 2])

            mars_trail.set_data(self.Mars.position[:frame, 0], self.Mars.position[:frame, 1])
            mars_trail.set_3d_properties(self.Mars.position[:frame, 2])

            saturn_trail.set_data(self.Saturn.position[:frame, 0], self.Saturn.position[:frame, 1])
            saturn_trail.set_3d_properties(self.Saturn.position[:frame, 2])

            uranus_trail.set_data(self.Uranus.position[:frame, 0], self.Uranus.position[:frame, 1])
            uranus_trail.set_3d_properties(self.Uranus.position[:frame, 2])

            neptune_trail.set_data(self.Neptune.position[:frame, 0], self.Neptune.position[:frame, 1])
            neptune_trail.set_3d_properties(self.Neptune.position[:frame, 2])

            return sun_plot, earth_plot, jupiter_plot, mars_plot, saturn_plot, earth_trail, jupiter_trail, mars_trail, saturn_trail, uranus_trail, neptune_trail

        ani = FuncAnimation(fig, update, frames=range(0, self.N,200), init_func=init, blit=False, interval=100)
        ani.save("Figures/solar_system_3D.gif", writer=PillowWriter(fps=30))
        
def main():
    sim = SolarSystem(N=int(365 * 100), dt=60 * 60 * 24)  # simulate for 100 years
    sim.simulate()
    sim.animate()

if __name__ == "__main__":
    main()
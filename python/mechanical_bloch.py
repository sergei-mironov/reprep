
import pennylane as qml
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
from numpy import ndarray
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Solution:
  t:ndarray
  x1:ndarray
  v1:ndarray
  x2:ndarray
  v2:ndarray


# Constants
m = 0.1  # mass in kg
l = 0.15  # length in m
k = 0.5  # spring constant in N/m
g = 9.81  # acceleration due to gravity in m/s^2

# Initial conditions
def coupled_oscillators(x1_0 = 0.01, x2_0 = 0.01, v1_0 = 0.0, v2_0 = 0.0):
  '''
  x1_0 = 0.01  # initial position of x1 in m
  x2_0 = 0.01  # initial position of x2 in m
  v1_0 = 0.0  # initial velocity of x1
  v2_0 = 0.0  # initial velocity of x2
  '''

  # System of equations
  def _ode(t, y) -> Solution:
    x1, v1, x2, v2 = y
    dx1dt = v1
    dv1dt = -(g / l) * x1 - (k / m) * (x1 - x2)
    dx2dt = v2
    dv2dt = -(g / l) * x2 - (k / m) * (x2 - x1)
    return [dx1dt, dv1dt, dx2dt, dv2dt]

  # Time span for the simulation
  t_span = (0, 10)  # simulate from t=0 to t=10 seconds
  t_eval = np.linspace(t_span[0], t_span[1], 1000)  # time points to evaluate at

  # Initial state vector
  y0 = [x1_0, v1_0, x2_0, v2_0]

  # Solve the system of ODEs
  solution = solve_ivp(_ode, t_span, y0, t_eval=t_eval)
  return Solution(solution.t, *[np.array(x) for x in solution.y])

def splot(name:str, sol:Solution)->None:
  # Plot the results on separate subplots
  plt.figure(figsize=(10, 8))

  # Plot for x1
  plt.subplot(211)
  plt.plot(sol.t, sol.x1, label=r'$x_1$', color='b')
  plt.ylabel('Displacement (m)')
  plt.title('Coupled Oscillator Simulation')
  plt.legend(loc='upper right')
  plt.grid(True)

  # Plot for x2
  plt.subplot(212, sharex=plt.gca())
  plt.plot(sol.t, sol.x2, label=r'$x_2$', color='r')
  plt.xlabel('Time (s)')
  plt.ylabel('Displacement (m)')
  plt.legend(loc='upper right')
  plt.grid(True)

  plt.tight_layout()
  plt.savefig(f"img/mechanical-bloch-f1-{name}.png")

# coupled_oscillators("1", 0.01, 0.01)
# coupled_oscillators("2", 0.01, -0.01)
# coupled_oscillators("3", 0.01, 0)

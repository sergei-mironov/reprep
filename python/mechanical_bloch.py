# import pennylane as qml
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
from numpy import ndarray
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from math import sqrt, cos

@dataclass
class PendulumProblem:
  m:float=0.1
  l:float=0.15
  k:float=0.5
  g:float=9.81

# {{{ Simulation

@dataclass
class Simulation:
  t:ndarray
  x1:ndarray
  v1:ndarray
  x2:ndarray
  v2:ndarray

# }}} Simulation


def coupled_pendulums(p:PendulumProblem, x10=0.01, x20=0.01, v10=0.0, v20=0.0) -> Simulation:
  """ Coupled pendulums without detuning. As described in "Waves and Oscillations. Prelude to
  Quantum Mechanincs".

  Arguments:
  - `x10`: initial position of x1 in m
  - `x20`: initial position of x2 in m
  - `v10`: initial velocity of x1
  - `v20`: initial velocity of x2
  """

  # Constants
  m, l, k, g = list(p.__dict__.values())

  # System of equations
  def _ode(t, y):
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
  y0 = [x10, v10, x20, v20]

  # Solve the system of ODEs
  solution = solve_ivp(_ode, t_span, y0, t_eval=t_eval)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])


# {{{ OscProblem

@dataclass
class OscProblem:
  m:float        = 0.9   # mass in kg
  k:float        = 5     # constant for main oscillating springs
  K:float        = 0.5   # constant for spring connecting two oscillators
  A:float        = 0.09  # FIXME: select appropriate
  sigma02:float  = (k + K) / m
  sigmac2:float  = K / m
  dsigma:float   = sigmac2 / sqrt(sigma02)
  wdrive:float   = dsigma
  delta:float    = dsigma - wdrive
  sigmaR:float   = sqrt(A**2 + delta**2)

# }}} OscProblem

def coupled_oscillators(p:OscProblem, xa0=0.01, xb0=0.01, va0=0.0, vb0=0.0)->Simulation:# {{{
  """ Coupled oscillators with periodic detuning, as described in the paper "The Classical Bloch
  Equations".
  """
  # Constants
  m, k, K, *_ = list(p.__dict__.values())

  nt = 1000 # number of time points to evaluate at

  # System of equations
  def _ode(t, state):
    xa, va, xb, vb = state
    dxadt = va
    dxbdt = vb
    dvadt = - xa * ((k + K) / m) + xb * (K / m)
    dvbdt = - xb * ((k + K) / m) + xa * (K / m)
    return [dxadt, dvadt, dxbdt, dvbdt]

  # Time span for the simulation
  t_span = (0, 10)  # simulate from t=0 to t=10 seconds
  t_eval = np.linspace(t_span[0], t_span[1], nt)

  # Initial state vector
  state0 = [xa0, va0, xb0, vb0]

  # Solve the system of ODEs
  solution = solve_ivp(_ode, t_span, state0, t_eval=t_eval)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])
# }}}

# {{{ coupled_detuned_oscillators

def coupled_detuned_oscillators(p:OscProblem, xa0=0.01, xb0=0.01, va0=0.0, vb0=0.0)->Simulation:
  """ Coupled oscillators with periodic detuning, as described in the paper "The Classical Bloch
  Equations". The function builds and sovles the system of ODEs corresponding to the problem.
  """
  m, k, K, *_ = list(p.__dict__.values())
  sigma02, A, wdrive = p.sigma02, p.A, p.wdrive
  sigma0 = sqrt(sigma02)
  nt = 1000
  t_span = (0, 90)
  t_eval = np.linspace(t_span[0], t_span[1], nt)
  state0 = [xa0, va0, xb0, vb0]

  def _ode(t, state):
    xa, va, xb, vb = state
    dk = -2.0 * sigma0 * m * A * cos(wdrive * t)
    dxadt = va
    dxbdt = vb
    dvadt = - xa * ((k + K) / m - dk / m) + xb * (K / m)
    dvbdt = - xb * ((k + K) / m + dk / m) + xa * (K / m)
    return [dxadt, dvadt, dxbdt, dvbdt]

  solution = solve_ivp(_ode, t_span, state0, t_eval=t_eval)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])

# }}} coupled_detuned_oscillators

def splot(name:str|None, sol:Simulation)->None:# {{{
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
  if name is None:
    plt.show()
  else:
    plt.savefig(f"img/mechanical-bloch-f1-{name}.png")
# }}}

def splotn(name:str|None, sol:Simulation)->None:
  # Plot the results on separate subplots
  plt.close()
  plt.figure(figsize=(10, 8))

  # Plot for x1
  plt.subplot(211)
  plt.plot(sol.t, sol.x1 + sol.x2, label=r'$x_+$', color='b')
  plt.ylabel('x+ (m)')
  plt.title('Coupled Oscillator Simulation')
  plt.legend(loc='upper right')
  plt.grid(True)

  # Plot for x2
  plt.subplot(212, sharex=plt.gca())
  plt.plot(sol.t, sol.x1 - sol.x2, label=r'$x_-$', color='r')
  plt.xlabel('Time (s)')
  plt.ylabel('x- (m)')
  plt.legend(loc='upper right')
  plt.grid(True)

  plt.tight_layout()
  if name is None:
    plt.show()
  else:
    plt.savefig(f"img/mechanical-bloch-f1-{name}.png")

# coupled_oscillators("1", 0.01, 0.01)
# coupled_oscillators("2", 0.01, -0.01)
# coupled_oscillators("3", 0.01, 0)

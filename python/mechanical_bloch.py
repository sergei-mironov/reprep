# import pennylane as qml
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
from numpy import ndarray
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from math import sqrt, cos, exp
from qutip import Bloch

# PendulumProblem {{{
@dataclass
class PendulumProblem:
  m:float=0.1
  l:float=0.15
  k:float=0.5
  g:float=9.81
# }}}

# {{{ Simulation

@dataclass
class Simulation:
  t  :ndarray               # Time ticks
  xa :ndarray               # Positions of the 1st bob
  va :ndarray               # Velocities of the 1st bob
  xb :ndarray               # Positions of the 2nd bob
  vb :ndarray               # Velocities of the 2nd bob
  xp :ndarray|None = None   # (+)-mode coordinates
  xm :ndarray|None = None   # (-)-mode coordinates
  xpA:ndarray|None = None   # (+)-mode amplitudes
  xmA:ndarray|None = None   # (-)-mode amplitudes
  def __post_init__(self):
    if self.xp is None:
      self.xp = self.xa + self.xb
    if self.xm is None:
      self.xm = self.xa - self.xb

# }}} Simulation

# {{{ Schedule

Time = float
DriveEnabled = list[tuple[Time,Time]]  # [(Start, Duration)]

class Schedule:
  """ Encodes the time input parameters. """
  tspan:tuple[Time,Time]   # Overall simulation time span, (Start, Stop).
  time:ndarray             # Time points within the simulation time span.
  drive:DriveEnabled       # Time segments when `wdrive` signal is enabled.

  def __init__(self, drive:DriveEnabled|None=None):
    nt = 1000
    self.tspan = (0, 90)
    self.time = np.linspace(self.tspan[0], self.tspan[1], nt)
    self.drive = drive if drive is not None else [(0.0,float('inf'))]

# }}} Schedule

def within(t:Time, sched:DriveEnabled) -> bool:# {{{
  for seg in sched:
    if seg[0] <= t and t < seg[0]+seg[1]:
      return True
  return False
# }}}

# class Initials {{{
@dataclass
class Initials:
  xa0 :float = 0.01
  va0 :float = 0.0
  xb0 :float = 0.01
  vb0 :float = 0.0
# }}}

def coupled_pendulums(p:PendulumProblem, s:Schedule, i:Initials) -> Simulation: # {{{
  """ Coupled pendulums without detuning. As described in "Waves and Oscillations. Prelude to
  Quantum Mechanincs".

  Arguments:
  - `xa0`: initial position of xa in m
  - `xb0`: initial position of xb in m
  - `va0`: initial velocity of xa
  - `vb0`: initial velocity of xb
  """

  # Constants
  m, l, k, g = list(p.__dict__.values())
  xa0, va0, xb0, vb0 = list(i.__dict__.values())

  # System of equations
  def _ode(t, y):
    xa, va, xb, vb = y
    dxadt = va
    dvadt = -(g / l) * xa - (k / m) * (xa - xb)
    dxbdt = vb
    dvbdt = -(g / l) * xb - (k / m) * (xb - xa)
    return [dxadt, dvadt, dxbdt, dvbdt]

  # Initial state vector
  y0 = [xa0, va0, xb0, vb0]

  # Solve the system of ODEs
  solution = solve_ivp(_ode, s.tspan, y0, t_eval=s.time)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])
# }}}

# {{{ OscProblem

@dataclass
class OscProblem:
  m:       float = 0.9   # mass in kg
  k:       float = 5     # constant for main oscillating springs
  K:       float = 0.5   # constant for spring connecting two oscillators
  A:       float = 0.09  # FIXME: select appropriate
  sigma02: float = (k + K) / m
  sigmac2: float = K / m
  dsigma:  float = sigmac2 / sqrt(sigma02)
  wdrive:  float = dsigma
  delta:   float = dsigma - wdrive
  sigmaR:  float = sqrt(A**2 + delta**2)

# }}} OscProblem

# {{{ coupled_oscillators

def coupled_oscillators(p:OscProblem, s:Schedule, i:Initials)->Simulation:
  """ Solve the Coupled oscillators without detuning, as described in the paper "The Classical Bloch
  Equations".
  """
  # Constants
  m, k, K, *_ = list(p.__dict__.values())
  xa0, va0, xb0, vb0 = list(i.__dict__.values())

  # System of equations
  def _ode(t, state):
    xa, va, xb, vb = state
    dxadt = va
    dxbdt = vb
    dvadt = - xa * ((k + K) / m) + xb * (K / m)
    dvbdt = - xb * ((k + K) / m) + xa * (K / m)
    return [dxadt, dvadt, dxbdt, dvbdt]

  # Initial state vector
  state0 = [xa0, va0, xb0, vb0]

  # Solve the system of ODEs
  solution = solve_ivp(_ode, s.tspan, state0, t_eval=s.time)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])

# }}}
# {{{ scheduled_coupled_detuned_oscillators

def scheduled_coupled_detuned_oscillators(p:OscProblem, s:Schedule, i:Initials)->Simulation:
  """ Solve the Coupled oscillators problem with detuning, as described in the "The Classical Bloch
  Equations" paper. Assume that detuning drive signal is enabled according to the schedule `s`.
  """
  xa0, va0, xb0, vb0 = list(i.__dict__.values())
  m, k, K, *_ = list(p.__dict__.values())
  sigma02, A = p.sigma02, p.A
  sigma0 = sqrt(sigma02)
  state0 = [xa0, va0, xb0, vb0]

  def _ode(t, state):
    xa, va, xb, vb = state
    dk = -2.0 * sigma0 * m * A * cos(p.wdrive * t) if within(t, s.drive) else 0.0
    dxadt = va
    dxbdt = vb
    dvadt = - xa * ((k + K) / m - dk / m) + xb * (K / m)
    dvbdt = - xb * ((k + K) / m + dk / m) + xa * (K / m)
    return [dxadt, dvadt, dxbdt, dvbdt]

  solution = solve_ivp(_ode, s.tspan, state0, t_eval=s.time)
  return Simulation(solution.t, *[np.array(x) for x in solution.y])

# }}} scheduled_coupled_detuned_oscillators
# {{{ coupled_detuned_oscillator_theoretic

def coupled_detuned_oscillator_theoretic(p:OscProblem, s:Schedule, i:Initials)->Simulation:
  assert s.drive and s.drive[0] == (0.0, float('inf'))
  t = s.time
  a0 = i.xa0 + i.xb0
  b0 = i.xa0 - i.xb0
  a_ = a0 * np.cos((p.A/2.0) * t) + 1j * b0 * np.sin((p.A/2.0) * t)
  b_ = b0 * np.cos((p.A/2.0) * t) + 1j * a0 * np.sin((p.A/2.0) * t)
  a = a_ * np.exp(-1j * (p.wdrive/2.0) * t)
  b = b_ * np.exp(+1j * (p.wdrive/2.0) * t)
  xp = a * np.exp(1j * np.sqrt(p.sigma02) * t)
  xm = b * np.exp(1j * np.sqrt(p.sigma02) * t)
  return Simulation(t, None, None, None, None,
                    xp=np.real(xp), xm=np.real(xm),
                    xpA=np.abs(xp), xmA=np.abs(xm))

# }}} coupled_detuned_oscillator_theoretic
# {{{ coupled_detuned_oscillators

def coupled_detuned_oscillators_vs_theoretic(p:OscProblem, i:Initials)->Simulation:
  """ Solve the coupled detuned oscillators problem numerically, set the normal mode amplitudes
  using theoretic solution (see the "Mechanical Bloch Equations" paper). Assume the detuning drive
  signal is always enabled.
  """
  s = Schedule()
  sim = scheduled_coupled_detuned_oscillators(p, s, i)
  sim_ref = coupled_detuned_oscillator_theoretic(p, s, i)
  return sim, sim_ref

# }}} coupled_detuned_oscillators

def splot(name:str|None, sol:Simulation)->None:# {{{
  """ Plot the simulation results, each oscialltor separately
  """

  # Plot the results on separate subplots
  plt.figure(figsize=(10, 8))

  # Plot for xa
  plt.subplot(211)
  plt.plot(sol.t, sol.xa, label=r'$x_1$', color='b')
  plt.ylabel('Displacement (m)')
  plt.title('Coupled Oscillator Simulation')
  plt.legend(loc='upper right')
  plt.grid(True)

  # Plot for xb
  plt.subplot(212, sharex=plt.gca())
  plt.plot(sol.t, sol.xb, label=r'$x_2$', color='r')
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
def splotn(name:str|None, sol:Simulation, sol_ref:Simulation|None=None)->str|None:# {{{
  """ Plot the coupled pendulim simulation normal mode results
  """
  plt.close("all")
  plt.figure(figsize=(10, 8))

  # Plot for xa
  plt.subplot(211)
  plt.plot(sol.t, sol.xp, label=r'Numeric $x_+$', color='b')
  if sol_ref is not None:
    plt.plot(sol_ref.t, sol_ref.xpA, label=r'Theoretic $|x_+|$', color='g')
  plt.ylabel('x+ (m)')
  plt.title('Coupled Oscillator Simulation')
  plt.legend(loc='upper right')
  plt.grid(True)

  # Plot for xb
  plt.subplot(212, sharex=plt.gca())
  plt.plot(sol.t, sol.xm, label=r'Numeric $x_-$', color='r')
  if sol_ref is not None:
    plt.plot(sol_ref.t, sol_ref.xmA, label=r'Theoretic $|x_-|$', color='g')
  plt.xlabel('Time (s)')
  plt.ylabel('x- (m)')
  plt.legend(loc='upper right')
  plt.grid(True)

  plt.tight_layout()
  if name is None:
    plt.show()
    return None
  else:
    f = f"img/mechanical-bloch-f1-{name}.png"
    plt.savefig(f)
    return f
# }}}
def splotb(name, p:OscProblem, s:Schedule, i:Initials):# {{{
  assert s.drive and s.drive[0] == (0.0, float('inf'))
  t = s.time
  a0 = i.xa0 + i.xb0
  b0 = i.xa0 - i.xb0

  # Calculate a_ and b_ with normalization
  a_ = a0 * np.cos((p.A / 2.0) * t) + 1j * b0 * np.sin((p.A / 2.0) * t)
  b_ = b0 * np.cos((p.A / 2.0) * t) + 1j * a0 * np.sin((p.A / 2.0) * t)
  norm_factor = np.sqrt(np.abs(a_)**2 + np.abs(b_)**2)

  # Normalizing a_ and b_
  a_ /= norm_factor
  b_ /= norm_factor

  # Calculate sx, sy, sz using the normalized a_ and b_
  sx = np.real(     a_ * np.conj(b_) +      np.conj(a_) * b_ )
  sy = np.real(1j * a_ * np.conj(b_) - 1j * np.conj(a_) * b_ )
  sz = np.real(     a_ * np.conj(a_) -      b_ * np.conj(b_) )

  # Stack sx, sy, sz into a single array with 3 rows
  points = np.stack((sx, sy, sz))
  print(points.shape)
  print(points[:,:4])


  # Set up the Bloch sphere
  b = Bloch(view=[-40,30])
  # Name the sphere poles as A and B
  b.zlabel = ['$A$', '$B$']
  b.add_points(points, 'l', 'r')
  b.show()
  if name is None:
    plt.show()
    return None
  else:
    f = f"img/mechanical-bloch-f1-{name}.png"
    plt.savefig(f)
    return f
# }}}

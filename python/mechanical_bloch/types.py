from math import sqrt, cos, exp
from dataclasses import dataclass

# class Initials: # {{{
@dataclass
class Initials:
  xa0 :float = 0.01
  va0 :float = 0.0
  xb0 :float = 0.01
  vb0 :float = 0.0
# }}}

# class OscProblem: # {{{
@dataclass
class OscProblem:
  m:       float = 0.9   # mass in kg
  k:       float = 5     # constant for main oscillating springs
  K:       float = 0.9   # constant for spring connecting two oscillators
  A:       float = 0.1   # amplitude of the spring constant oscillation
  sigma02: float = (k + K) / m
  sigmac2: float = K / m
  dsigma:  float = sigmac2 / sqrt(sigma02)
  wdrive:  float = dsigma
  delta:   float = dsigma - wdrive
  sigmaR:  float = sqrt(A**2 + delta**2)
# }}}

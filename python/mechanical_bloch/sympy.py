import sympy as sp
from math import pi
from .types import Initials


def bless_complex(v):
  if not isinstance(v, sp.Basic) or not v.is_complex:
    raise ValueError(f"{v} is not a complex number.")
  return complex(v)

# Define the constants
m = 0.9   # mass in kg
k = 5     # constant for main oscillating springs
K = 0.9   # constant for spring connecting two oscillators
A = 0.1   # amplitude of the spring constant oscillation

# Calculate derived constants
sigma02 = (k + K) / m
sigmac2 = K / m
dsigma = sigmac2 / sp.sqrt(sigma02)
wdrive = dsigma
delta = dsigma - wdrive
sigmaR = sp.sqrt(A**2 + delta**2)

# Define a0 and b0 as complex constants
a0 = sp.symbols('a0', complex=True)
b0 = sp.symbols('b0', complex=True)

def _calc(t):
  # FIXME: probably one need to put `a0_ = a0 + b0` and `b0_ = a0 - b0`.
  # Calculate expressions for a_ and b_
  a_ = a0 * sp.cos((A/2.0) * t) + sp.I * b0 * sp.sin((A/2.0) * t)
  b_ = b0 * sp.cos((A/2.0) * t) + sp.I * a0 * sp.sin((A/2.0) * t)

  # Compute a, b
  a = a_ * sp.exp(-sp.I * (wdrive/2.0) * t)
  b = b_ * sp.exp(sp.I * (wdrive/2.0) * t)

  # Compute xp, xm
  xp = a * sp.exp(sp.I * sp.sqrt(sigma02) * t)
  xm = b * sp.exp(sp.I * sp.sqrt(sigma02) * t)

  # Compute xa, xb
  xa = (xp + xm) / 2
  xb = (xp - xm) / 2
  return xa, xb

# Represent time with a symbolic variable
t = sp.symbols('t', real=True)

# Calculate xa and xb using the calc function
xa, xb = _calc(t)

def xa_(t_, a0_, b0_):
  return bless_complex(xa.subs({
    t: t_, a0: a0_, b0: b0_
  }))

def xb_(t_, a0_, b0_):
  return bless_complex(xb.subs({
    t: t_, a0: a0_, b0: b0_
  }))


# Calculate the expressions for the time-derivatives va and vb
va = sp.diff(xa, t)
vb = sp.diff(xb, t)


def va_(t_, a0_, b0_):
  return bless_complex(va.subs({
    t: t_, a0: a0_, b0: b0_
  }))

def vb_(t_, a0_, b0_):
  return bless_complex(vb.subs({
    t: t_, a0: a0_, b0: b0_
  }))



print('va',va)
print('vb',vb)

# Define the initial conditions xa0, xb0, va0, vb0 as symbols
xa0, xb0, va0, vb0 = sp.symbols('xa0 xb0 va0 vb0', real=True)

# Solve for a0 and b0 using the initial conditions
initial_eqs = (
  sp.Eq(sp.re(xa.subs(t, 0)), xa0),  # xa(t=0) = xa0
  sp.Eq(sp.re(xb.subs(t, 0)), xb0),  # xb(t=0) = xb0
  sp.Eq(sp.re(va.subs(t, 0)), va0),  # va(t=0) = va0
  sp.Eq(sp.re(vb.subs(t, 0)), vb0)   # vb(t=0) = vb0
)

solutions = sp.solve(initial_eqs, (a0, b0))
assert len(solutions) == 1
solution = solutions[0]
print(solution[a0])
print(solution[b0])

def a0_(initials: Initials) -> complex:
  # Substitute values from the Initials dataclass
  a0_value = solution[a0].subs({
    xa0: initials.xa0,
    va0: initials.va0,
    xb0: initials.xb0,
    vb0: initials.vb0
  })
  return bless_complex(a0_value)

def b0_(initials: Initials) -> complex:
  # Substitute values from the Initials dataclass
  b0_value = solution[b0].subs({
    xa0: initials.xa0,
    va0: initials.va0,
    xb0: initials.xb0,
    vb0: initials.vb0
  })
  return bless_complex(b0_value)


def test(d=pi/4):
  i = Initials()
  t1, t2 = pi/2 - d, pi/2 + d
  print('xa', [xa_(10*t, a0_(i), b0_(i)) for t in [t1, t2]])
  print('va', [va_(10*t, a0_(i), b0_(i)) for t in [t1, t2]])
  print('xb', [xb_(10*t, a0_(i), b0_(i)) for t in [t1, t2]])
  print('vb', [vb_(10*t, a0_(i), b0_(i)) for t in [t1, t2]])

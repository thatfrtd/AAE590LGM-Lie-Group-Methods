import sympy as sym
import casadi as cas

# Sympy
R = sym. MatrixSymbol('R', 3, 3)
eta = sym.MatrixSymbol('eta', 3, 1)
Omega = sym.MatrixSymbol('Omega', 3, 1)

phi = R.T @ eta

print(phi.as_explicit().jacobian(eta))

# CasADi
R = cas.SX.sym('R', 3, 3)
eta = cas.SX.sym('eta', 3, 1)
Omega = cas.SX.sym('Omega', 3, 1)

phi = R.T @ eta

print(cas.jacobian(phi, eta))

# Frechet derivative

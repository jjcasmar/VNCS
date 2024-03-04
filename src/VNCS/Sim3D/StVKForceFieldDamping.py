import sympy as sym
from sympy import diff

F00 = sym.Symbol('F00')
F01 = sym.Symbol('F01')
F10 = sym.Symbol('F10')
F11 = sym.Symbol('F11')
F20 = sym.Symbol('F20')
F21 = sym.Symbol('F21')

dFdt00 = sym.Symbol('F00')
dFdt01 = sym.Symbol('F01')
dFdt10 = sym.Symbol('F10')
dFdt11 = sym.Symbol('F11')
dFdt20 = sym.Symbol('F20')
dFdt21 = sym.Symbol('F21')

F = sym.Matrix([[F00, F01], [F10, F11], [F20, F21]])
dFdt = sym.Matrix([[dFdt00, dFdt01], [dFdt10, dFdt11], [dFdt20, dFdt21]])

e = 1.0 / 2.0 * (F.T @ dFdt + dFdt.T @ F)

e00 = e[0, 0]
e11 = e[1, 1]
e10 = e[1, 0]

E = e00 * e00 + e11 * e11 + e10 * e10

diff_E = sym.Matrix([
    [sym.diff(E, dFdt00), sym.diff(E, dFdt01)],
    [sym.diff(E, dFdt10), sym.diff(E, dFdt11)],
    [sym.diff(E, dFdt20), sym.diff(E, dFdt21)],
])


# diff2_E_X = [
#    [diff(diff_E[0, 0], F00), diff(diff_E[0, 0], F01), diff(diff_E[0, 0], F10), diff(
#        diff_E[0, 0], F11), diff(diff_E[0, 0], F20), diff(diff_E[0, 0], F21)],
#    [diff(diff_E[0, 1], F00), diff(diff_E[0, 1], F01), diff(diff_E[0, 1], F10), diff(
#        diff_E[0, 1], F11), diff(diff_E[0, 1], F20), diff(diff_E[0, 1], F21)],
#    [diff(diff_E[1, 0], F00), diff(diff_E[1, 0], F01), diff(diff_E[1, 0], F10), diff(
#        diff_E[1, 0], F11), diff(diff_E[1, 0], F20), diff(diff_E[1, 0], F21)],
#    [diff(diff_E[1, 1], F00), diff(diff_E[1, 1], F01), diff(diff_E[1, 1], F10), diff(
#        diff_E[1, 1], F11), diff(diff_E[1, 1], F20), diff(diff_E[1, 1], F21)],
#    [diff(diff_E[2, 0], F00), diff(diff_E[2, 0], F01), diff(diff_E[2, 0], F10), diff(
#        diff_E[2, 0], F11), diff(diff_E[2, 0], F20), diff(diff_E[2, 0], F21)],
#    [diff(diff_E[2, 1], F00), diff(diff_E[2, 1], F01), diff(diff_E[2, 1], F10), diff(
#        diff_E[2, 1], F11), diff(diff_E[2, 1], F20), diff(diff_E[2, 1], F21)]
# ]
#
# diff2_E_V = [
#    [diff(diff_E[0, 0], dFdt00), diff(diff_E[0, 0], dFdt01), diff(diff_E[0, 0], dFdt10), diff(
#        diff_E[0, 0], dFdt11), diff(diff_E[0, 0], dFdt20), diff(diff_E[0, 0], dFdt21)],
#    [diff(diff_E[0, 1], dFdt00), diff(diff_E[0, 1], dFdt01), diff(diff_E[0, 1], dFdt10), diff(
#        diff_E[0, 1], dFdt11), diff(diff_E[0, 1], dFdt20), diff(diff_E[0, 1], dFdt21)],
#    [diff(diff_E[1, 0], dFdt00), diff(diff_E[1, 0], dFdt01), diff(diff_E[1, 0], dFdt10), diff(
#        diff_E[1, 0], dFdt11), diff(diff_E[1, 0], dFdt20), diff(diff_E[1, 0], dFdt21)],
#    [diff(diff_E[1, 1], dFdt00), diff(diff_E[1, 1], dFdt01), diff(diff_E[1, 1], dFdt10), diff(
#        diff_E[1, 1], dFdt11), diff(diff_E[1, 1], dFdt20), diff(diff_E[1, 1], dFdt21)],
#    [diff(diff_E[2, 0], dFdt00), diff(diff_E[2, 0], dFdt01), diff(diff_E[2, 0], dFdt10), diff(
#        diff_E[2, 0], dFdt11), diff(diff_E[2, 0], dFdt20), diff(diff_E[2, 0], dFdt21)],
#    [diff(diff_E[2, 1], dFdt00), diff(diff_E[2, 1], dFdt01), diff(diff_E[2, 1], dFdt10), diff(
#        diff_E[2, 1], dFdt11), diff(diff_E[2, 1], dFdt20), diff(diff_E[2, 1], dFdt21)]
# ]

print(diff_E)

# print(diff2_E_X)
# print(diff2_E_V)

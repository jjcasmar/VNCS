syms F00 F01 F10 F11 F20 F21 real;
syms dFdt00 dFdt01 dFdt10 dFdt11 dFdt20 dFdt21 real;

F = [F00 F01; F10 F11; F20 F21];
dFdt = [dFdt00 dFdt01; dFdt10 dFdt11; dFdt20 dFdt21];

e = 1 / 2 * (F' *dFdt + dFdt' * F);

E = e(1, 1) * e(1, 1) + e(2, 2) * e(2, 2) + e(2, 1) * e(2, 1);

dE = [
  diff(E, dFdt00) diff(E, dFdt01); diff(E, dFdt10) diff(E, dFdt11);
  diff(E, dFdt20), diff(E, dFdt21)
];

dE_X = [
    diff(dE(1,1), F00), diff(dE(1,1), F01), diff(dE(1,1), F10), diff(dE(1,1), F11), diff(dE(1,1), F20), diff(dE(1,1), F21);
    diff(dE(1,2), F00), diff(dE(1,2), F01), diff(dE(1,2), F10), diff(dE(1,2), F11), diff(dE(1,2), F20), diff(dE(1,2), F21);
    diff(dE(2,1), F00), diff(dE(2,1), F01), diff(dE(2,1), F10), diff(dE(2,1), F11), diff(dE(2,1), F20), diff(dE(2,1), F21);
    diff(dE(2,2), F00), diff(dE(2,2), F01), diff(dE(2,2), F10), diff(dE(2,2), F11), diff(dE(2,2), F20), diff(dE(2,2), F21);
    diff(dE(3,1), F00), diff(dE(3,1), F01), diff(dE(3,1), F10), diff(dE(3,1), F11), diff(dE(3,1), F20), diff(dE(3,1), F21);
    diff(dE(3,2), F00), diff(dE(3,2), F01), diff(dE(3,2), F10), diff(dE(3,2), F11), diff(dE(3,2), F20), diff(dE(3,2), F21);
];

dE_V= [
    diff(dE(1,1), dFdt00), diff(dE(1,1), dFdt01), diff(dE(1,1), dFdt10), diff(dE(1,1), dFdt11), diff(dE(1,1), dFdt20), diff(dE(1,1), dFdt21);
    diff(dE(1,2), dFdt00), diff(dE(1,2), dFdt01), diff(dE(1,2), dFdt10), diff(dE(1,2), dFdt11), diff(dE(1,2), dFdt20), diff(dE(1,2), dFdt21);
    diff(dE(2,1), dFdt00), diff(dE(2,1), dFdt01), diff(dE(2,1), dFdt10), diff(dE(2,1), dFdt11), diff(dE(2,1), dFdt20), diff(dE(2,1), dFdt21);
    diff(dE(2,2), dFdt00), diff(dE(2,2), dFdt01), diff(dE(2,2), dFdt10), diff(dE(2,2), dFdt11), diff(dE(2,2), dFdt20), diff(dE(2,2), dFdt21);
    diff(dE(3,1), dFdt00), diff(dE(3,1), dFdt01), diff(dE(3,1), dFdt10), diff(dE(3,1), dFdt11), diff(dE(3,1), dFdt20), diff(dE(3,1), dFdt21);
    diff(dE(3,2), dFdt00), diff(dE(3,2), dFdt01), diff(dE(3,2), dFdt10), diff(dE(3,2), dFdt11), diff(dE(3,2), dFdt20), diff(dE(3,2), dFdt21);
];
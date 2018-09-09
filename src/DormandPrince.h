#pragma once

#include "KerrSchild.h"

// dp_derivs - Berechnung eines RK-Schrittes
void dp_derivs(Kerrschild &H, double* R, double* R_Out);

// rkdp - Dormand-Price-Verfahren
void dopri(Kerrschild &Hole, double* y, double* dydx, double &dt, double &T, double err_tol);

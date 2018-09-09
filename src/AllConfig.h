#pragma once

#include <cmath>

static constexpr double PIov2 = 1.57079632679489661923;

//Anzahl Pixel in Breite und Hoehe
static constexpr int PIXEL_WIDTH = 1600;
static constexpr int PIXEL_HEIGHT = 900;

// Einstellungen Integration
static constexpr double dopri_error = 0.00000001;
static constexpr double T_MAX = 4.0;

//Daten zu schwarzem Loch
/* Drehung der Rotationsachse des Schwarzen Lochs
* yz_tilt: dreht bei Blick entgegen x-Achse entgegen dem Uhrzeigersinn
* xy_tilt: dreht bei Blick entgegen z-Achse entgegen dem Uhrzeigersinn
*/
static constexpr int yz_tilt = 45;
static constexpr int xy_tilt = -45;

static constexpr double M = 0.02;
static constexpr double a = 0.99 * M;
static const double EnteredHole = M + sqrt(M*M - a*a);
static const double EnteredHoleSq = EnteredHole * EnteredHole;
static constexpr double xHole = 0;
static constexpr double yHole = 0;
static constexpr double zHole = 0;

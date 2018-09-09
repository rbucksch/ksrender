#pragma once

#include "Camera.h"
#include "Wall.h"

// TraceRay - Berechnet alle Endpunkte des Ray-Tracings
void traceRays(Camera &Cam, Wall Bounds[]);

// ObstacleHit - ueberprueft, ob eine Grenze getroffen wurde
int ObstacleImpact(double* now, double* last, Wall Bounds[]);

#pragma once

#include "AllConfig.h"

#include <cmath>


class Kerrschild
{
private:
	static double rs;
public:
	Kerrschild();
	~Kerrschild();

	void tiltDown(double* input);
	void tiltUp(double* input);

	double gammaScale(double* V);

	void updateGamma(double* V);
	void update(double* V);

	double alpha;

	double dLNalphadx;
	double dLNalphady;
	double dLNalphadz;

	double dalphadx;
	double dalphady;
	double dalphadz;

	double beta_x;
	double beta_y;
	double beta_z;

	double dbeta_xdx;
	double dbeta_xdy;
	double dbeta_xdz;
	double dbeta_ydx;
	double dbeta_ydy;
	double dbeta_ydz;
	double dbeta_zdx;
	double dbeta_zdy;
	double dbeta_zdz;

	double ga_xx;
	double ga_xy;
	double ga_xz;
	double ga_yy;
	double ga_yz;
	double ga_zz;

	double ga_inv_xx;
	double ga_inv_xy;
	double ga_inv_xz;
	double ga_inv_yy;
	double ga_inv_yz;
	double ga_inv_zz;

	double K_xx;
	double K_xy;
	double K_xz;
	double K_yy;
	double K_yz;
	double K_zz;

	double g_x_xx;
	double g_x_xy;
	double g_x_xz;
	double g_x_yy;
	double g_x_yz;
	double g_x_zz;

	double g_y_xx;
	double g_y_xy;
	double g_y_xz;
	double g_y_yy;
	double g_y_yz;
	double g_y_zz;

	double g_z_xx;
	double g_z_xy;
	double g_z_xz;
	double g_z_yy;
	double g_z_yz;
	double g_z_zz;
};

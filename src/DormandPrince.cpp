#include "AllConfig.h"
#include "DormandPrince.h"
#include "KerrSchild.h"

#include <iostream>
#include <cmath>
#include <iomanip>

void dp_derivs(Kerrschild &H, double* R, double* R_Out)
{
	double u_x = R[5];
	double u_y = R[6];
	double u_z = R[7];

	H.update(R);

	double K_x = H.K_xx*u_x + H.K_xy*u_y + H.K_xz*u_z;
	double K_y = H.K_xy*u_x + H.K_yy*u_y + H.K_yz*u_z;
	double K_z = H.K_xz*u_x + H.K_yz*u_y + H.K_zz*u_z;

	double Kx_x = H.ga_inv_xx*H.K_xx + H.ga_inv_xy*H.K_xy + H.ga_inv_xz*H.K_xz;
	double Kx_y = H.ga_inv_xx*H.K_xy + H.ga_inv_xy*H.K_yy + H.ga_inv_xz*H.K_yz;
	double Kx_z = H.ga_inv_xx*H.K_xz + H.ga_inv_xy*H.K_yz + H.ga_inv_xz*H.K_zz;
	double Ky_x = H.ga_inv_xy*H.K_xx + H.ga_inv_yy*H.K_xy + H.ga_inv_yz*H.K_xz;
	double Ky_y = H.ga_inv_xy*H.K_xy + H.ga_inv_yy*H.K_yy + H.ga_inv_yz*H.K_yz;
	double Ky_z = H.ga_inv_xy*H.K_xz + H.ga_inv_yy*H.K_yz + H.ga_inv_yz*H.K_zz;
	double Kz_x = H.ga_inv_xz*H.K_xx + H.ga_inv_yz*H.K_xy + H.ga_inv_zz*H.K_xz;
	double Kz_y = H.ga_inv_xz*H.K_xy + H.ga_inv_yz*H.K_yy + H.ga_inv_zz*H.K_yz;
	double Kz_z = H.ga_inv_xz*H.K_xz + H.ga_inv_yz*H.K_yz + H.ga_inv_zz*H.K_zz;

	double G_xx = H.g_x_xx*u_x + H.g_x_xy*u_y + H.g_x_xz*u_z;
	double G_xy = H.g_x_xy*u_x + H.g_x_yy*u_y + H.g_x_yz*u_z;
	double G_xz = H.g_x_xz*u_x + H.g_x_yz*u_y + H.g_x_zz*u_z;
	double G_yx = H.g_y_xx*u_x + H.g_y_xy*u_y + H.g_y_xz*u_z;
	double G_yy = H.g_y_xy*u_x + H.g_y_yy*u_y + H.g_y_yz*u_z;
	double G_yz = H.g_y_xz*u_x + H.g_y_yz*u_y + H.g_y_zz*u_z;
	double G_zx = H.g_z_xx*u_x + H.g_z_xy*u_y + H.g_z_xz*u_z;
	double G_zy = H.g_z_xy*u_x + H.g_z_yy*u_y + H.g_z_yz*u_z;
	double G_zz = H.g_z_xz*u_x + H.g_z_yz*u_y + H.g_z_zz*u_z;

	R_Out[0] = -1;
	R_Out[1] = -H.alpha*u_x + H.beta_x;
	R_Out[2] = -H.alpha*u_y + H.beta_y;
	R_Out[3] = -H.alpha*u_z + H.beta_z;
	R_Out[4] = 0;
	R_Out[5] = -H.alpha * (u_x*(u_x*(H.dLNalphadx - K_x) + 2.0*Kx_x - G_xx)
						 + u_y*(u_x*(H.dLNalphady - K_y) + 2.0*Kx_y - G_xy)
						 + u_z*(u_x*(H.dLNalphadz - K_z) + 2.0*Kx_z - G_xz))
						 + H.ga_inv_xx*H.dalphadx
						 + H.ga_inv_xy*H.dalphady
						 + H.ga_inv_xz*H.dalphadz
						 + H.dbeta_xdx*u_x
						 + H.dbeta_xdy*u_y
						 + H.dbeta_xdz*u_z;

	R_Out[6] = -H.alpha * (u_x*(u_y*(H.dLNalphadx - K_x) + 2.0*Ky_x - G_yx)
						 + u_y*(u_y*(H.dLNalphady - K_y) + 2.0*Ky_y - G_yy)
						 + u_z*(u_y*(H.dLNalphadz - K_z) + 2.0*Ky_z - G_yz))
						 + H.ga_inv_xy*H.dalphadx
						 + H.ga_inv_yy*H.dalphady
						 + H.ga_inv_yz*H.dalphadz
						 + H.dbeta_ydx*u_x
						 + H.dbeta_ydy*u_y
						 + H.dbeta_ydz*u_z;

	R_Out[7] = -H.alpha * (u_x*(u_z*(H.dLNalphadx - K_x) + 2.0*Kz_x - G_zx)
						 + u_y*(u_z*(H.dLNalphady - K_y) + 2.0*Kz_y - G_zy)
						 + u_z*(u_z*(H.dLNalphadz - K_z) + 2.0*Kz_z - G_zz))
						 + H.ga_inv_xz*H.dalphadx
						 + H.ga_inv_yz*H.dalphady
						 + H.ga_inv_zz*H.dalphadz
						 + H.dbeta_zdx*u_x
						 + H.dbeta_zdy*u_y
						 + H.dbeta_zdz*u_z;
}

void dopri(Kerrschild &H, double* y, double* dydx, double &dt, double &T, double err_tol)
{
	static const double c21 = 0.2;				static const double c31 = 0.075;
	static const double c32 = 0.225;			static const double c41 = 44.0/45.0;
	static const double c42 = -56.0/15.0;		static const double c43 = 32.0/9.0;
	static const double c51 = 19372.0/6561.0;	static const double c52 = -25360.0/2187.0;
	static const double c53 = 64448.0/6561.0;	static const double c54 = -212.0/729.0;
	static const double c61 = 9017.0/3168.0;	static const double c62 = -355.0/33.0;
	static const double c63 = 46732.0/5247.0;	static const double c64 = 49.0/176.0;
	static const double c65 = -5103.0/18656.0;	static const double c71 = 35.0/384.0;
	static const double c73 = 500.0/1113.0;		static const double c74 = 125.0/192.0;
	static const double c75 = -2187.0/6784.0;	static const double c76 = 11.0/84.0;
	
	static const double b1 = 5179.0/57600.0;	static const double b3 = 7571.0/16695.0;
	static const double b4 = 393.0/640.0;		static const double b5 = -92097.0/339200.0;
	static const double b6 = 187.0/2100.0;		static const double b7 = 0.025;

	static const double e1 = 71.0/57600.0;		static const double e3 = -71.0/16695.0;
	static const double e4 = 71.0/1920.0;		static const double e5 = -17253.0/339200.0;
	static const double e6 = 22.0/525.0;		static const double e7 = -0.025;
	
	double k2[8];
	double k3[8];
	double k4[8];
	double k5[8];
	double k6[8];
	double dydx_new[8];
	double w[8];

	double p[4];
	double o[8];

	//dydx already calculated; equivalent to k1
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*c21*dydx[i];
	}

	dp_derivs(H, w, k2);
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*(c31*dydx[i] + c32*k2[i]);
	}

	dp_derivs(H, w, k3);
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*(c41*dydx[i] + c42*k2[i] + c43*k3[i]);
	}

	dp_derivs(H, w, k4);
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*(c51*dydx[i] + c52*k2[i] + c53*k3[i] + c54*k4[i]);
	}

	dp_derivs(H, w, k5);
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*(c61*dydx[i] + c62*k2[i] + c63*k3[i] + c64*k4[i] + c65*k5[i]);
	}
	
	dp_derivs(H, w, k6);

	//Endwerte berechnen
	for (int i = 0; i<8; i++) {
		w[i] = y[i] + dt*(c71*dydx[i] + c73*k3[i] + c74*k4[i] + c75*k5[i] + c76*k6[i]);
	}

	//k7 des aktuellen & k1 des naechsten Schrittes berechnen
	dp_derivs(H, w, dydx_new);

	//Fehler abschaetzen
	for (int i = 1; i<4; i++) {
		p[i] = y[i] + dt*(b1*dydx[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*dydx_new[i]);
	}
	for (int i = 1; i<8; i++) {
		o[i] = dt*(e1*dydx[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*dydx_new[i]);
	}

	H.updateGamma(p);
	double errx = H.gammaScale(o + 1);
	double errv = H.gammaScale(o + 5);

	double err1 = (errv > errx ? errv : errx);

	if (err1 < err_tol) {
		for (int i = 0; i<8; i++) {
			y[i] = w[i];
			dydx[i] = dydx_new[i];
		}
		T += dt;
	}

	double delta = 0.9 * pow(err_tol / err1, 0.2);

	if (delta >= 5) {
		dt = dt * 5;
	}
	else if (delta <= 0.1) {
		dt = dt * 0.1;
	}
	else {
		dt = delta * dt;
	}

	if (dt < 0.00000000001) {
		dt = 0.00000000001;
	}
	else if (dt > 0.1) {
		dt = 0.1;
	}

}

#include "KerrSchild.h"
#include "AllConfig.h"

#include <cmath>

static const double dxy = xy_tilt * PIov2 / 90.0;
static const double dyz = yz_tilt * PIov2 / 90.0;

static const double trig[8] = { cos(dxy), sin(dxy), cos(dyz)*sin(dxy), cos(dyz)*cos(dxy), sin(dyz), sin(dyz)*sin(dxy), sin(dyz)*cos(dxy), cos(dyz) };

Kerrschild::Kerrschild()
{
}

Kerrschild::~Kerrschild()
{
}

void Kerrschild::tiltDown(double* input)
{
	double x = input[1];
	double y = input[2];
	double z = input[3];

	double vx = input[5];
	double vy = input[6];
	double vz = input[7];

	input[1] = trig[0] * x + trig[1] * y;
	input[2] = -trig[2] * x + trig[3] * y + trig[4] * z;
	input[3] = trig[5] * x - trig[6] * y + trig[7] * z;

	input[5] = trig[0] * vx + trig[1] * vy;
	input[6] = -trig[2] * vx + trig[3] * vy + trig[4] * vz;
	input[7] = trig[5] * vx - trig[6] * vy + trig[7] * vz;
}

void Kerrschild::tiltUp(double* input)
{
	double x = input[1];
	double y = input[2];
	double z = input[3];

	double vx = input[5];
	double vy = input[6];
	double vz = input[7];

	input[1] = trig[0] * x - trig[2] * y + trig[5] * z;
	input[2] = trig[1] * x + trig[3] * y - trig[6] * z;
	input[3] = trig[4] * y + trig[7] * z;

	input[5] = trig[0] * vx - trig[2] * vy + trig[5] * vz;
	input[6] = trig[1] * vx + trig[3] * vy - trig[6] * vz;
	input[7] = trig[4] * vy + trig[7] * vz;
}

double Kerrschild::gammaScale(double* V)
{
	double d = sqrt(ga_xx*V[0]*V[0]	+ 2.0*ga_xy*V[0]*V[1] + 2.0*ga_xz*V[0]*V[2]
									+	  ga_yy*V[1]*V[1] + 2.0*ga_yz*V[1]*V[2]
														  +		ga_zz*V[2]*V[2] );
	return d;
}

void Kerrschild::updateGamma(double* V)
{
	double x = V[1];
	double y = V[2];
	double z = V[3];

	double a2 = a * a;
	double z2 = z * z;

	double R2 = x*x + y*y + z2 - a2;
	double r2 = 0.5*(R2 + sqrt(R2*R2 + 4.0*a2*z2));
	double r = sqrt(r2);

	double ra = 1.0 / (r2 + a2);

	double f = r*r2 / (rs*(r2*r2 + a2*z2));

	double lx = ra*(x*r + a*y);
	double ly = ra*(y*r - a*x);
	double lz = z / r;

	ga_xx = f * lx*lx + 1.0;
	ga_xy = f * lx*ly;
	ga_xz = f * lx*lz;
	ga_yy = f * ly*ly + 1.0;
	ga_yz = f * ly*lz;
	ga_zz = f * lz*lz + 1.0;
}

void Kerrschild::update(double* V)
{
	double x = V[1];
	double y = V[2];
	double z = V[3];

	double a2 = a * a;
	double z2 = z * z;

	double R2 = x*x + y*y + z2 - a2;
	double r2 = 0.5*(R2 + sqrt(R2*R2 + 4.0*a2*z2));
	double r = sqrt(r2);

	double ra = 1.0 / (r2 + a2);
	double rInv = 1.0 / r;
	double r2Inv = 1.0 / r2;

	double f = r*r2 / (rs*(r2*r2 + a2*z2));

	double lx = ra * (x*r + a*y);
	double ly = ra * (y*r - a*x);
	double lz = z*rInv;
	
	double frs = f*rs;
	double drdx = frs * x;
	double drdy = frs * y;
	double drdz = frs * z * (1.0 + a2*r2Inv);

	double dfconst = f*(3.0*rInv - 4.0*frs);
	double dfdx = dfconst*drdx;
	double dfdy = dfconst*drdy;
	double dfdz = dfconst*drdz - 2.0*a*a*f*frs*lz*r2Inv;

	double dlxconst = ra*(x - 2.0*lx*r);
	double dlxdx = dlxconst*drdx + r*ra;
	double dlxdy = dlxconst*drdy + a*ra;
	double dlxdz = dlxconst*drdz;

	double dlyconst = ra*(y - 2.0*ly*r);
	double dlydx = dlyconst*drdx - a*ra;
	double dlydy = dlyconst*drdy + r*ra;
	double dlydz = dlyconst*drdz;


	double dlzdx = -drdx*lz*rInv;
	double dlzdy = -drdy*lz*rInv;
	double dlzdz = (1.0 - drdz*lz)*rInv;

	double alpha2 = 1.0 / (f + 1.0);

	double b_faktor = f * alpha2;
	double db_faktor = alpha2 * alpha2;
	double g_faktor = 0.5 * alpha2;

	alpha = sqrt(alpha2);

	double bfOVr = b_faktor*rInv;
	double dLNconst1 = bfOVr*(2.0*f*r*rs - 1.5);
	double dLNconst2 = dLNconst1 + bfOVr;

	dLNalphadx = dLNconst1*drdx;
	dLNalphady = dLNconst1*drdy;
	dLNalphadz = dLNconst2*drdz - bfOVr*frs*z;

	beta_x = b_faktor * lx;
	beta_y = b_faktor * ly;
	beta_z = b_faktor * lz;

	dalphadx = alpha * dLNalphadx;
	dalphady = alpha * dLNalphady;
	dalphadz = alpha * dLNalphadz;

	double ff1 = f*(f + 1.0);
	dbeta_xdx = (dfdx*lx + dlxdx*ff1)*db_faktor;
	dbeta_xdy = (dfdy*lx + dlxdy*ff1)*db_faktor;
	dbeta_xdz = (dfdz*lx + dlxdz*ff1)*db_faktor;
	dbeta_ydx = (dfdx*ly + dlydx*ff1)*db_faktor;
	dbeta_ydy = (dfdy*ly + dlydy*ff1)*db_faktor;
	dbeta_ydz = (dfdz*ly + dlydz*ff1)*db_faktor;
	dbeta_zdx = (dfdx*lz + dlzdx*ff1)*db_faktor;
	dbeta_zdy = (dfdy*lz + dlzdy*ff1)*db_faktor;
	dbeta_zdz = (dfdz*lz + dlzdz*ff1)*db_faktor;

	ga_xx = f * lx*lx + 1.0;
	ga_xy = f * lx*ly;
	ga_xz = f * lx*lz;
	ga_yy = f * ly*ly + 1.0;
	ga_yz = f * ly*lz;
	ga_zz = f * lz*lz + 1.0;

	ga_inv_xx = alpha2 * (f + 2.0 - ga_xx);
	ga_inv_xy = -alpha2 * ga_xy;
	ga_inv_xz = -alpha2 * ga_xz;
	ga_inv_yy = alpha2 * (f + 2.0 - ga_yy);
	ga_inv_yz = -alpha2 * ga_yz;
	ga_inv_zz = alpha2 * (f + 2.0 - ga_zz);


	K_xx = 0.5*alpha * (lx*(dfdx*(ga_xx + 1.0) + dfdz * ga_xz + dfdy * ga_xy) + 2.0*f*(ga_xz*(dlxdz - dlzdx) + ga_xy * (dlxdy - dlydx) + dlxdx));

	K_xy = 0.5*alpha * (ga_xx*dfdx*ly + ga_yy * dfdy*lx + ga_xy * dfdz*lz + f * ((ga_xx - ga_yy)*(dlydx - dlxdy) + dlydx + dlxdy + (dlxdz - dlzdx)*ga_yz + (dlydz - dlzdy)*ga_xz));

	K_xz = 0.5*alpha * (ga_xx*dfdx*lz + ga_zz * dfdz*lx + ga_xz * dfdy*ly + f * ((ga_xx - ga_zz)*(dlzdx - dlxdz) + dlzdx + dlxdz + (dlxdy - dlydx)*ga_yz + (dlzdy - dlydz)*ga_xy));

	K_yy = 0.5*alpha * (ly*(dfdy*(ga_yy + 1.0) + dfdz * ga_yz + dfdx * ga_xy) + 2.0*f*(ga_yz*(dlydz - dlzdy) + ga_xy * (dlydx - dlxdy) + dlydy));

	K_yz = 0.5*alpha * (ga_yy*dfdy*lz + ga_zz * dfdz*ly + ga_yz * dfdx*lx + f * ((ga_yy - ga_zz)*(dlzdy - dlydz) + dlzdy + dlydz + (dlydx - dlxdy)*ga_xz + (dlzdx - dlxdz)*ga_xy));

	K_zz = 0.5*alpha * (lz*(dfdz*(ga_zz + 1.0) + dfdy * ga_yz + dfdx * ga_xz) + 2.0*f*(ga_yz*(dlzdy - dlydz) + ga_xz * (dlzdx - dlxdz) + dlzdz));


	g_x_xx = lx * (2.0*(f*dlxdx + ga_xy * (dlxdy - dlydx) + ga_xz * (dlxdz - dlzdx)) + lx * (dfdx*(ga_xx - f) + dfdy * ga_xy + dfdz * ga_xz))*g_faktor;

	g_x_xy = lx * (-dfdx * ly*(ga_yy + ga_zz - 2.0) + lx * (dfdy*ga_yy + dfdz * ga_yz) + f * (ga_yz*(dlxdz - dlzdx) + ga_xz * (dlydz - dlzdy) + (dlxdy - dlydx)*(ga_zz + 2.0*ga_yy - 1.0) + 2.0*dlydx))*g_faktor;

	g_x_xz = lx * (-dfdx * lz*(ga_yy + ga_zz - 2.0) + lx * (dfdy*ga_yz + dfdz * ga_zz) + f * (ga_yz*(dlxdy - dlydx) + ga_xy * (dlzdy - dlydz) + (dlxdz - dlzdx)*(ga_yy + 2.0*ga_zz - 1.0) + 2.0*dlzdx))*g_faktor;

	g_x_yy = (lx*(2.0*f*(ga_yz*(dlydz - dlzdy) + dlydy) + ly * (dfdz*ga_yz + dfdy * (ga_yy + 1.0))) + ly * (ga_xx - f - 2.0)*(dfdx*ly + 2.0*f*(dlydx - dlxdy)))*g_faktor;

	g_x_yz = (lx*(f*((dlzdy - dlydz)*(ga_yy - ga_zz) + dlzdy + dlydz) + dfdz * ly*ga_zz + dfdy * lz*ga_yy) + (ga_xx - f - 2)*(dfdx*ly*lz + f * ((dlydx - dlxdy)*lz + (dlzdx - dlxdz)*ly)))*g_faktor;

	g_x_zz = (lx*(2.0*f*(ga_yz*(dlzdy - dlydz) + dlzdz) + lz * (dfdy*ga_yz + dfdz * (ga_zz + 1.0))) + lz * (ga_xx - f - 2.0)*(dfdx*lz + 2.0*f*(dlzdx - dlxdz)))*g_faktor;

	g_y_xx = (ly*(2.0*f*(ga_xz*(dlxdz - dlzdx) + dlxdx) + lx * (dfdz*ga_xz + dfdx * (ga_xx + 1.0))) + lx * (ga_yy - f - 2.0)*(dfdy*lx + 2.0*f*(dlxdy - dlydx)))*g_faktor;

	g_y_xy = ly * (-dfdy * lx*(ga_xx + ga_zz - 2.0) + ly * (dfdx*ga_xx + dfdz * ga_xz) + f * (ga_xz*(dlydz - dlzdy) + ga_yz * (dlxdz - dlzdx) + (dlydx - dlxdy)*(2.0*ga_xx + ga_zz - 1.0) + 2.0*dlxdy))*g_faktor;

	g_y_xz = (ly*(f*((dlzdx - dlxdz)*(ga_xx - ga_zz) + dlzdx + dlxdz) + dfdz * lx*ga_zz + dfdx * lz*ga_xx) + (ga_yy - f - 2)*(dfdy*lx*lz + f * ((dlxdy - dlydx)*lz + (dlzdy - dlydz)*lx)))*g_faktor;

	g_y_yy = ly * (2.0*(f*dlydy + ga_xy * (dlydx - dlxdy) + ga_yz * (dlydz - dlzdy)) + ly * (dfdy*(ga_yy - f) + dfdx * ga_xy + dfdz * ga_yz))*g_faktor;

	g_y_yz = ly * (-dfdy * lz*(ga_xx + ga_zz - 2.0) + ly * (dfdx*ga_xz + dfdz * ga_zz) + f * (ga_xz*(dlydx - dlxdy) + ga_xy * (dlzdx - dlxdz) + (dlydz - dlzdy)*(ga_xx + 2.0*ga_zz - 1.0) + 2.0*dlzdy))*g_faktor;

	g_y_zz = (ly*(2.0*f*(ga_xz*(dlzdx - dlxdz) + dlzdz) + lz * (dfdx*ga_xz + dfdz * (ga_zz + 1.0))) + lz * (ga_yy - f - 2.0)*(dfdy*lz + 2.0*f*(dlzdy - dlydz)))*g_faktor;

	g_z_xx = (lz*(2.0*f*(ga_xy*(dlxdy - dlydx) + dlxdx) + lx * (dfdy*ga_xy + dfdx * (ga_xx + 1.0))) + lx * (ga_zz - f - 2.0)*(dfdz*lx + 2.0*f*(dlxdz - dlzdx)))*g_faktor;

	g_z_xy = (lz*(f*((dlydx - dlxdy)*(ga_xx - ga_yy) + dlydx + dlxdy) + dfdy * lx*ga_yy + dfdx * ly*ga_xx) + (ga_zz - f - 2)*(dfdz*lx*ly + f * ((dlxdz - dlzdx)*ly + (dlydz - dlzdy)*lx)))*g_faktor;


	g_z_xz = lz * (-dfdz * lx*(ga_xx + ga_yy - 2.0) + lz * (dfdx*ga_xx + dfdy * ga_xy) + f * (ga_xy*(dlzdy - dlydz) + ga_yz * (dlxdy - dlydx) - (dlxdz - dlzdx)*(2.0*ga_xx + ga_yy - 1.0) + 2.0*dlxdz))*g_faktor;

	g_z_yy = (lz*(2.0*f*(ga_xy*(dlydx - dlxdy) + dlydy) + ly * (dfdx*ga_xy + dfdy * (ga_yy + 1.0))) + ly * (ga_zz - f - 2.0)*(dfdz*ly + 2.0*f*(dlydz - dlzdy)))*g_faktor;

	g_z_yz = lz * (lz * (dfdx * ga_xy + dfdy * ga_yy) - dfdz * ly*(ga_xx + ga_yy - 2.0) - f * (ga_xy*(dlxdz - dlzdx) + ga_xz * (dlxdy - dlydx) + (dlydz - dlzdy)*(ga_xx + 2.0*ga_yy - 1.0) - 2.0*dlydz))*g_faktor;

	g_z_zz = lz * (2.0*(f*dlzdz + ga_yz * (dlzdy - dlydz) + ga_xz * (dlzdx - dlxdz)) + lz * (dfdz*(ga_zz - f) + dfdy * ga_yz + dfdx * ga_xz))*g_faktor;
}

double Kerrschild::rs = 0.5 / M;
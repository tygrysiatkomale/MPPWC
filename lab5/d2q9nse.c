#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vis.h"

#define F 9
#define cs_sq (1./3.)
#define wC (4./9.)
#define wS (1./9.)
#define wL (1./36.)

#define AR(x,y,f) (f + F*(x) + F*LX*(y))
#define ARU(x,y) (x + LX*(y))

enum {_C, _N, _S, _W, _E, _NE, _NW, _SW, _SE};
enum {_FLUID, _SOLID, _INLET, _OUTLETX};

const int LX = (2*128);
const int LY = (1*128);

double cx[F] = {0, 0, 0, -1, 1, 1,-1,-1, 1};
double cy[F] = {0, 1, -1, 0, 0, 1, 1,-1,-1};
double w[F] = {wC, wS, wS, wS, wS, wL, wL, wL, wL};

double *cells;
double *temp_cells;

double tau;
double *uLB, *vLB, *rhoLB;
int *map;
double rho_IN, ux_IN, uy_IN;

double fEq(int i, double rho, double u, double v) {
    double cu = (cx[i]*u + cy[i]*v)/cs_sq;
    return w[i] * rho * (1 + cu + 0.5*cu*cu - 0.5*(u*u +v*v)/cs_sq);
}

void collideAndStream(double *c, double *tc) {
 int x, y, i, xp, yp, xm, ym, mp;
 double fstar[F];
 double rho,ux,uy;
 
    for (y = 0; y < LY; y++) {
        for (x = 0; x < LX; x++) {
            mp = map[ARU(x,y)];
            // BC
            if (mp == _SOLID) {
                fstar[_C] = c[AR(x,y,_C)];
                fstar[_N] = c[AR(x,y,_S)];
                fstar[_S] = c[AR(x,y,_N)];
                fstar[_W] = c[AR(x,y,_E)];
                fstar[_E] = c[AR(x,y,_W)];
                fstar[_NE] = c[AR(x,y,_SW)];
                fstar[_SW] = c[AR(x,y,_NE)];
                fstar[_NW] = c[AR(x,y,_SE)];
                fstar[_SE] = c[AR(x,y,_NW)];

                rhoLB[ARU(x,y)] = 1.;
                uLB[ARU(x,y)] = 0.;
                vLB[ARU(x,y)] = 0.;
            } else if ( mp == _INLET ) {
                for( i = 0; i < F; i++ )
                    fstar[i] = fEq(i, rho_IN, ux_IN, uy_IN);

                rhoLB[ARU(x,y)] = rho_IN;
                uLB[ARU(x,y)] = ux_IN;
                vLB[ARU(x,y)] = uy_IN;                    
            } else if ( mp == _OUTLETX ) {
                double rho_xm = rhoLB[ARU(x-1,y)]; // 1.
                double ux_xm = uLB[ARU(x-1,y)];
                double uy_xm = vLB[ARU(x-1,y)];
                for( i = 0; i < F; i++ ) {
                    fstar[i] = fEq(i, rho_xm, ux_xm, uy_xm);
                }
                rhoLB[ARU(x,y)] = rho_xm;
                uLB[ARU(x,y)] = ux_xm;
                vLB[ARU(x,y)] = uy_xm;                    
            } else {// Compute macroscopic variables
                rho = 0; ux = 0; uy = 0;
                for (i = 0; i < F; i++) {
                    rho += c[AR(x,y,i)];
                    ux += cx[i]*c[AR(x,y,i)];
                    uy += cy[i]*c[AR(x,y,i)];
                }
                ux /= rho;
                uy /= rho;

                rhoLB[ARU(x,y)] = rho;
                uLB[ARU(x,y)] = ux;
                vLB[ARU(x,y)] = uy;
                // Collision step
                // MRT
                double m_rho,m_eps,m_e,m_jx,m_jy,m_qx,m_qy,m_pxx,m_pxy;
                double m_rhoe,m_epse,m_ee,m_jxe,m_jye,m_qxe,m_qye,m_pxxe,m_pxye;
                double f_C,f_N,f_S,f_W,f_E,f_NE,f_NW,f_SW,f_SE;
                f_C = c[AR(x,y,_C)];
                f_N = c[AR(x,y,_N)];
                f_S = c[AR(x,y,_S)];
                f_W = c[AR(x,y,_W)];
                f_E = c[AR(x,y,_E)];
                f_NE = c[AR(x,y,_NE)];
                f_NW = c[AR(x,y,_NW)];
                f_SW = c[AR(x,y,_SW)];
                f_SE = c[AR(x,y,_SE)];
                // m = Mf momenty
                m_rho = rho;
                m_jx = rho*ux;
                m_jy = rho*uy;
                m_e = -(4*f_C + f_N + f_S + f_W + f_E) + 2*(f_NE + f_NW + f_SW + f_SE);
                m_eps = 4*f_C -2*(f_N + f_S + f_W + f_E) + f_NE + f_NW + f_SW + f_SE;
                m_qx = 2*(f_W-f_E) + f_NE + f_SE - f_NW - f_SW;
                m_qy = 2*(f_S-f_N) + f_NW - f_SW + f_NE - f_SE;
                m_pxx = f_E + f_W - f_N - f_S;
                m_pxy = f_NE - f_NW - f_SE + f_SW;
                // rownowagi momentow
                m_rhoe = m_rho;
                m_jxe = m_jx;
                m_jye = m_jy;
                m_ee = -2*rho + 3*(m_jx*m_jx + m_jy*m_jy);
                m_epse = rho - 3*(m_jx*m_jx + m_jy*m_jy);
                m_qxe = -m_jx;
                m_qye = -m_jy;
                m_pxxe = m_jx*m_jx-m_jy*m_jy;
                m_pxye = m_jx*m_jy;
                // czestosci relaksacji
                double om_e, om_eps, om_q, om_nu;
                om_e = 1.63;
                om_eps = 1.54;
                om_q = 1.9;
                om_nu = 1./tau;
                // kolizja momentow
                m_e = om_e*(m_ee - m_e);
                m_eps = om_eps*(m_epse - m_eps);
                m_qx = om_q*(m_qxe - m_qx);
                m_qy = om_q*(m_qye - m_qy);
                m_pxx = om_nu*(m_pxxe - m_pxx);
                m_pxy = om_nu*(m_pxye - m_pxy);
                m_rho = m_jx = m_jy = 0;

                fstar[_C] = f_C + 1./9.*(m_eps - m_e);
                fstar[_E] = f_E + 1./36.*(-m_e - 2*m_eps - 6*m_qx + 9*m_pxx);
                fstar[_N] = f_N + 1./36.*(-m_e - 2*m_eps - 6*m_qy - 9*m_pxx);
                fstar[_W] = f_W + 1./36.*(-m_e - 2*m_eps + 6*m_qx + 9*m_pxx);
                fstar[_S] = f_S + 1./36.*(-m_e - 2*m_eps + 6*m_qy - 9*m_pxx);
                fstar[_NE] = f_NE + 1./36.*(2*m_e + m_eps + 3*m_qx + 3*m_qy + 9*m_pxy);
                fstar[_NW] = f_NW + 1./36.*(2*m_e + m_eps - 3*m_qx + 3*m_qy - 9*m_pxy);
                fstar[_SW] = f_SW + 1./36.*(2*m_e + m_eps - 3*m_qx - 3*m_qy + 9*m_pxy);
                fstar[_SE] = f_SE + 1./36.*(2*m_e + m_eps + 3*m_qx - 3*m_qy - 9*m_pxy);

                // for (i = 0; i < F; i++) {
                //     fstar[i] = (1-1./tau)*c[AR(x,y,i)] + 
                //             fEq(i, rho, ux, uy)/tau;
                //     // if (fstar[i] < 0 || fstar[i] > 2*w[i]) {
                //     //     fstar[i] = w[i]; // Prevent negative values
                //     // }
                // }
            }
            // Streaming step & BC
            xp = (x + 1) % LX;
            xm = (x - 1 + LX) % LX;
            yp = (y + 1) % LY;
            ym = (y - 1 + LY) % LY;

            tc[AR(xp,yp,_NE)] = fstar[_NE];
            tc[AR(xm,yp,_NW)] = fstar[_NW];
            tc[AR(xm,ym,_SW)] = fstar[_SW];
            tc[AR(xp,ym,_SE)] = fstar[_SE];
            tc[AR(xp,y,_E)] = fstar[_E];
            tc[AR(xm,y,_W)] = fstar[_W];
            tc[AR(x,yp,_N)] = fstar[_N];
            tc[AR(x,ym,_S)] = fstar[_S];
            tc[AR(x,y,_C)] = fstar[_C];
        }
    }
}

void 
setInitialConditions(double rhoi, double ui, double vi, double *c) {
 int x, y, i;
    for (y = 0; y < LY; y++) {
        for (x = 0; x < LX; x++) {

            if ( y == 0 || y == LY-1 ) map[ARU(x,y)] = _SOLID;
            else if (x == 0) map[ARU(x,y)] = _INLET;
            else if (x == LX-1) map[ARU(x,y)] = _OUTLETX;
            else map[ARU(x,y)] = _FLUID;

            // Set velocity to double sheer layer pattern
            uLB[ARU(x,y)] = ui;// * sin(2*M_PI*x/LX) + cos(2*M_PI*y/LY);
            vLB[ARU(x,y)] = vi;// * sin(2*M_PI*x/LX) + cos(2*M_PI*y/LY);
            rhoLB[ARU(x,y)] = rhoi;
            for (i = 0; i < F; i++) {
                c[AR(x,y,i)] = fEq(i, rhoi, uLB[ARU(x,y)], vLB[ARU(x,y)]);
            }
        }
    }
    
 for ( y = LY/2 - 8; y < LY/2 + 8; y++) 
    for (x = 4 + LX/4 - 8; x < 4 + LX/4 + 8; x++)
        map[ARU(x,y)] = _SOLID;  

 for ( y = 8 + LY/2 - 8; y < 8 +LY/2 + 8; y++) 
    for (x = 30 + LX/4 - 8; x < 30 + LX/4 + 8; x++)
        map[ARU(x,y)] = _SOLID;  


}

void dumpStateVTK(char *fname) {
 FILE *fp;
 int x,y;

 fp = fopen(fname, "w");
 fprintf(fp,"# vtk DataFile Version 2.0\n");
 fprintf(fp,"2D-ADE data file \n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET RECTILINEAR_GRID\n");
 fprintf(fp,"DIMENSIONS %d %d %d\n", LX, LY, 1);
 fprintf(fp,"X_COORDINATES %d int\n", LX);
 for (x = 0; x < LX; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Y_COORDINATES %d int\n", LY);
 for (x = 0; x < LY; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Z_COORDINATES 1 int\n");
 fprintf(fp, "0\n");
 fprintf(fp,"POINT_DATA %d \n", LX*LY);
 fprintf(fp,"SCALARS pressure double 1\n");
 fprintf(fp,"LOOKUP_TABLE default\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++) {
    double rho = rhoLB[ARU(x,y)];
    fprintf(fp, "%e \n", rho);
  }
 fprintf(fp,"SCALARS velocity_mag double 1\n");
 fprintf(fp,"LOOKUP_TABLE default\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++) {
    double velmag = hypot(uLB[ARU(x,y)],vLB[ARU(x,y)]);
    fprintf(fp, "%e \n", velmag);
  }fprintf(fp,"\n");
  fprintf(fp,"VECTORS velocity double\n");
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++) {
    fprintf(fp, "%e %e 0.\n", uLB[ARU(x,y)], vLB[ARU(x,y)]);
  }fprintf(fp,"\n");
 fclose(fp);
}


int main(int argc, char *argv[]) {
int iter = 0;
int maxIter = 1000000;

double Ma, Re, N;

Ma = 0.1;
Re = 100000;
N = LY;

double rhoi = 1;
double ui = 0; 
double vi = 0;

rho_IN = 1;
ux_IN = Ma*sqrt(cs_sq);
uy_IN = 0;

double nuLB = ux_IN*N/Re;

cells = (double*)malloc(LX*LY*F*sizeof(double));
temp_cells = (double*)malloc(LX*LY*F*sizeof(double));

uLB = (double*)malloc(LX*LY*sizeof(double));
vLB = (double*)malloc(LX*LY*sizeof(double));
rhoLB = (double*)malloc(LX*LY*sizeof(double));
map = (int*)malloc(LX*LY*sizeof(int));

tau = nuLB/cs_sq + 0.5;

printf("tau = %f\n", tau);
printf("nuLB = %f\n", nuLB);
printf("ux_IN = %f\n", ux_IN);

setInitialConditions(rhoi, ui, vi, cells);
// visualization
vis_t rho_vis, vel_vis;

vis_init(&rho_vis, LX, LY, 4, "rho");
vis_init(&vel_vis, LX, LY, 4, "|u|");

vis_draw_rho(&rho_vis, rhoLB,uLB, vLB, LX, LY, F);
vis_draw_speed(&vel_vis, rhoLB, uLB, vLB, LX, LY, F);

dumpStateVTK("output_000.vtk");
do {
    if (iter % 2)
        collideAndStream(temp_cells, cells);
    else
        collideAndStream(cells, temp_cells);
    iter++;

    vis_draw_rho(&rho_vis, rhoLB, uLB, vLB, LX, LY, F);
    vis_draw_speed(&vel_vis, rhoLB, uLB, vLB, LX, LY, F);

} while (iter < maxIter);
dumpStateVTK("output_final.vtk");

free(cells);
free(temp_cells);

vis_close(&rho_vis);
vis_close(&vel_vis);

return 0;
}
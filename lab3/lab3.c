#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define F 9
#define cs_sq (1./3.)
#define wC (4./9.)
#define wS (1./9.)
#define wL (1./36.)
#define M_PI 3.14159265358979323846

#define AR(x, y, f) (f + F*(x) + F*LX*(y))
#define ARU(x, y) (x + LX*(y))

enum {_C, _N, _S, _W, _E, _NE, _NW, _SW, _SE};

const int LX = 128;
const int LY = 128;

double cx[F] = {0, 0, 0, -1, 1, 1, -1, -1, 1};
double cy[F] = {0, 1, -1, 0, 0, 1, 1, -1, -1};
double w[F] = {wC, wS, wS, wS, wS, wL, wL, wL, wL};

double *cells;
double *temp_cells;

double tau;
double *uLB, *vLB;
double rho_IN, ux_IN, uy_IN;

double fEq(int i, double phi, double u, double v){
    double cu = (cx[i]*u + cy[i]*v)/cs_sq;
    return w[i] * phi * (1 + cu + 0.5*cu*cu - 0.5*(u*u + v*v)/cs_sq); 
}

void collideAndStream(double *c, double *tc) {
    int x, y, i, xp, yp, xm, ym;
    double fstar[F];
    double rho, ux, uy;
    double tmp;

    for (y =0 ; y < LY ; y++) {
        for (x = 0 ; x < LX ; x++) {
            //BC
            if ( x == 0 || y == 0 || x == LX - 1) {
                tmp = c[AR(x, y, _W)];
                c[AR(x, y, _W)] = c[AR(x, y, _E)];
                c[AR(x, y, _E)] = tmp;
                tmp = c[AR(x, y, _N)];
                c[AR(x, y, _N)] = c[AR(x, y, _S)];
                c[AR(x, y, _S)] = tmp;
                tmp = c[AR(x, y, _NW)];
                c[AR(x, y, _NW)] = c[AR(x, y, _SE)];
                c[AR(x, y, _SE)] = tmp;
                tmp = c[AR(x, y, _NE)];
                c[AR(x, y, _NE)] = c[AR(x, y, _SW)];
                c[AR(x, y, _SW)] = tmp;
                for (i = 0 ; i < F ; i++) {
                        fstar[i] = c[AR(x,y,i)];
                    }
            } else {
                // BC top
                if ( y == LY -1 ) {
                    for (i = 0 ; i < F ; i++) {
                        c[AR(x, y, i)] = fEq(i, rho_IN, ux_IN, uy_IN);
                    }
                }
                //Compute macroscopic variables
                rho = 0;
                ux = 0;
                uy = 0;
                for (i = 0 ; i < F ; i++) {
                    rho += c[AR(x, y, i)];
                    ux += cx[i] * c[AR(x, y, i)];
                    uy += cy[i] * c[AR(x, y, i)];
                }
                ux /= rho;
                uy /= rho;
                uLB[ARU(x,y)] = ux;
                vLB[ARU(x,y)] = uy;
                // Collision step
                for (i = 0 ; i < F ; i++) {
                    fstar[i] = (1-1./tau)*c[AR(x,y,i)] + 
                        fEq(i, rho, ux, uy)/tau;
                        if (fstar[i] < 0 || fstar[i] > w[i]) {
                            fstar[i] = w[i];
                        }
                }
            }
            // Streaming step & BC
            xp = (x + 1) % LX;
            xm = (x - 1 + LX) % LX;
            yp = (y + 1) % LY;
            ym = (y - 1 + LY) % LY;

            tc[AR(xp, yp, _NE)] = fstar[_NE];
            tc[AR(xm, ym, _SW)] = fstar[_SW];
            tc[AR(xm, yp, _NW)] = fstar[_NW];
            tc[AR(xp, ym, _SE)] = fstar[_SE];
            tc[AR(x, yp, _N)] = fstar[_N];
            tc[AR(x, ym, _S)] = fstar[_S];
            tc[AR(xm, y, _W)] = fstar[_W];
            tc[AR(xp, y, _E)] = fstar[_E];
            tc[AR(x, y, _C)] = fstar[_C];
        }
    }
}

void setInitialConditions(double rhoi, double ui, double vi, double *c) {
    int x, y, i;
    for (y = 0 ; y < LY ; y++) {
        for (x = 0 ; x < LX ; x++) {
            // Set the velocity to Taylor-Green vortex
            uLB[ARU(x,y)] = ui; //* sin(2*M_PI*x/LX) + cos(2*M_PI*y/LY);
            vLB[ARU(x,y)] = vi; //* sin(2*M_PI*x/LX) + cos(2*M_PI*y/LY);
            for (i = 0 ; i < F ; i++) {
                c[AR(x, y, i)] = fEq(i, rhoi, uLB[ARU(x,y)], vLB[ARU(x,y)]);
            }
        }
    }
    
    // for ( y = LY/2 ; y < LY/2 + 8 ; y++) {
    //     for (x = LX/2 ; x < LX/2 + 8 ; x++) {
    //         for (i = 0 ; i < F ; i++) {
    //             c[AR(x, y, i)] = fEq(i, 1.0, ui, vi);
    //         }
    //     }
    // }
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
    fprintf(fp,"SCALARS temperature double 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (y = 0; y < LY; y++)
    for (x = 0; x < LX; x++) {
        double velmag = hypot(uLB[ARU(x,y)], vLB[ARU(x,y)]);
        fprintf(fp, "%e \n", velmag);
    }
    fprintf(fp,"\n");
    fclose(fp);
}

int main(int argc, char *argv[]) {
    int iter = 0;
    int maxIter = 10;

    double Ma, Re, N;

    Ma = 0.1;
    Re = 1000000;
    N = LX;

    double rhoi = 1;
    double ui = 0;
    double vi = 0;

    rho_IN = 1;
    ux_IN = Ma*sqrt(cs_sq);
    uy_IN = 0;

    double nuLB = ux_IN * N / Re;

    cells = (double*)malloc(LX*LY*F*sizeof(double));
    temp_cells = (double*)malloc(LX*LY*F*sizeof(double));

    tau = nuLB / cs_sq + 0.5;

    printf("tau: %f\n", tau);
    printf("nuLB: %f\n", nuLB);
    printf("ux_IN: %f\n", ux_IN);

    uLB = (double*)malloc(LX*LY*sizeof(double));
    vLB = (double*)malloc(LX*LY*sizeof(double));
    
    setInitialConditions(rhoi, ui, vi, cells);
    dumpStateVTK("output_000.vtk");

    do {
        if (iter % 2) {
            collideAndStream(temp_cells, cells);
        } else {
            collideAndStream(cells, temp_cells);
        }
        iter++;
    } while (iter < maxIter);

    dumpStateVTK("output_final.vtk");

    free(cells);
    free(temp_cells);

    return 0;
}
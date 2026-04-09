#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define F 5
#define cs_sq (1./3.)
#define wC (1./3.)
#define wNC (1./6.)

#define AR(x, y, f) (f + F*(x) + F*LX*(y))

enum {_C, _N, _S, _W, _E};

const int LX = 128;
const int LY = 128;

double cx[F] = {0, 0, 0, -1, 1};
double cy[F] = {0, 1, -1, 0, 0};
double w[F] = {wC, wNC, wNC, wNC, wNC};

double *cells;
double *temp_cells;

double tau, uLB, vLB;

double fEq(int i, double phi, double u, double v){
    return w[i] * phi * (1 + (cx[i]*u + cy[i]*v)/cs_sq); 
}

void collideAndStream(double *c, double *tc) {
    int x, y, i, xp, yp, xm, ym;
    double fstar[F];
    double phi;

    for (y =0 ; y < LY ; y++) {
        for (x = 0 ; x < LX ; x++) {
            phi = 0;
            for (i = 0 ; i < F ; i++) {
                phi += c[AR(x, y, i)];;
            }
            // Collision step
            for (i = 0 ; i < F ; i++) {
                fstar[i] = c[AR(x,y,i)] + 1./tau*(fEq(i, phi, uLB, vLB) - c[AR(x,y,i)]);
            }
            // Streaming step & BC
            xp = (x + 1) % LX;
            xm = (x - 1 + LX) % LX;
            yp = (y + 1) % LY;
            ym = (y - 1 + LY) % LY;

            tc[AR(xp, y, _E)] = fstar[_E];
            tc[AR(xm, y, _W)] = fstar[_W];
            tc[AR(x, yp, _N)] = fstar[_N];
            tc[AR(x, ym, _S)] = fstar[_S];
            tc[AR(x, y, _C)] = fstar[_C];
        }
    }
}

void setInitialConditions(double phii, double ui, double vi, double *c) {
    int x, y, i;
    for (y = 0 ; y < LY ; y++) {
        for (x = 0 ; x < LX ; x++) {
            for (i = 0 ; i < F ; i++) {
                c[AR(x, y, i)] = fEq(i, phii, ui, vi);
            }
        }
    }

    for (i = 0 ; i < F ; i++) {
        c[AR(LX/2, LY/2, i)] = fEq(i, 1.0, ui, vi);
    }
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
        double rho = 0;
        int f;
        for (f = 0; f < F; f++) rho += cells[AR(x,y,f)];
        fprintf(fp, "%e \n", rho);
    }
    fprintf(fp,"\n");
    fclose(fp);
}

int main(int argc, char *argv[]) {
    int iter = 0;
    int maxIter = 1000;

    double phii = 0;
    double ui = 0;
    double vi = 0;

    double D = 0.01;

    cells = (double*)malloc(LX*LY*F*sizeof(double));
    temp_cells = (double*)malloc(LX*LY*F*sizeof(double));

    tau = D/cs_sq + 0.5;
    
    setInitialConditions(phii, ui, vi, cells);
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
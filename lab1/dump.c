// ParaView
// VisIt

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


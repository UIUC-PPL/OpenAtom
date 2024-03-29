#include "matrix2file.h"
#include "charm++.h"//< Just for CkAbort!!
#include <assert.h>

void dumpMatrixUber(int u, const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{

  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);

}

void dumpMatrixUberIter(int u, int iter, const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{

  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,".i.%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, iter, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);

}


void dumpMatrix(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);
}


void dumpMatrixUber(int u, const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g %.12g\n",x+i,y+j,matrix[i*ydim+j].re, matrix[i*ydim+j].im);
  fclose(loutfile);
}

void dumpMatrixUberIter(int u, int iter, const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,".i.%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, iter, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g %.12g\n",x+i,y+j,matrix[i*ydim+j].re, matrix[i*ydim+j].im);
  fclose(loutfile);
}



void dumpMatrix(const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g %.12g\n",x+i,y+j,matrix[i*ydim+j].re, matrix[i*ydim+j].im);
  fclose(loutfile);
}



void loadMatrix(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j])));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}

void loadMatrixUber(int u, const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j])));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}


void loadMatrixUberIter(int u, int iter,const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,".i.%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, iter, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  //  CkPrintf("[%d] loading %s\n",CkMyPe(),filename);
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j])));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}



void loadMatrix(const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf %lf\n",&junk1,&junk2,&(matrix[i*ydim+j].re),&(matrix[i*ydim+j].im)));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}

void loadMatrixUber(int u, const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf %lf\n",&junk1,&junk2,&(matrix[i*ydim+j].re),&(matrix[i*ydim+j].im)));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}

void loadMatrixUberIter(int u, int iter, const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,".i.%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, iter, u, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  //  CkPrintf("[%d] loading %s\n",CkMyPe(), filename);
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=0;j<ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf %lf\n",&junk1,&junk2,&(matrix[i*ydim+j].re),&(matrix[i*ydim+j].im)));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}


//! NOTE: this uses the evil piny convention
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=1;j<=ydim;j++)
        assert(fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i][j])));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}


//! NOTE: this uses the evil piny convention
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %d\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  memset(fmt, 0 , 1000);
  memset(filename, 0 , 1000);
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
  {
    int junk1,junk2;
    for(int i=0;i<xdim;i++)
      for(int j=1;j<=ydim;j++)
        assert(fscanf(loutfile,"%d %d %d\n",&junk1,&junk2,&(matrix[i][j])));
    fclose(loutfile);
  }
  else
  {
    CkAbort(filename);
  }
}


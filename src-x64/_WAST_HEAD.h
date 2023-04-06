#ifndef WAST_H_INCLUDED
#define WAST_H_INCLUDED

#define GAMMA_M 5
#define MPI 3.1415926
#define MPI1 0.1591549   //1.0/(2*MPI)
#define MPI2 0.3989423   //1.0/sqrt(2*MPI)
#define SQRT2 0.7071068  //1.0/sqrt(2.0)
#define LOG2PI 1.837877  //log(2pi)
#define MEPS 1e-10
#define ITMAX 100
#define EPS 1.0e-8
#define FPMIN 1.0e-30
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define IDEX(a, b) (((a)<(b))?1:0)

void AbyB(double *outVector, double *A, double *v, int n, int p, int q);

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q);

double exp_approx(double x, int n);

void Standarize(double *y, double *std, int n, int p, double *x, int flag);

void Std(double *std, double *x, int n, int p, int flag);

void sortN(int *ind0, double *x, int n, int dd);

void SortQ(double *s, int l, int r);

int LowTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

void SampleQuantile(double *qr, int m, double *z, int n, double *q);

int MatrixInvSymmetric(double *a,int n);

int CholeskyDecomL(double *L, int n, double *a);

int CholeskyDecomU(double *U, int n, double *a);

double incbeta(double x, double a, double b);

void Kernelh(double *x, double *kern, int n, double h0, int type);

void EstQR2(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps);

void EstLogisticR(double *beta, double *residual, double *x, double *y, int n, int p, int maxstep, double eps);

void EstLinearR(double *beta, double *residual, const double *x, const double *y, int n, int p);

double EstQuantileR1(double *y, double *residual, double tau, const int n);

void EstQuantileR(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps);

void Est_probitR(int *y, double *x, double *alpha, double *residual, int n, int p, int maxIter, double tol);

#endif // WAST_H_INCLUDED

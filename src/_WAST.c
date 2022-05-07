#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

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

void AbyB(double *outVector, double *A, double *v, int n, int p, int q){
	int i,j,k;
	double tmp;
	for (i=0;i<n;i++){
		for(k=0;k<q;k++){
			tmp = 0;
			for(j=0;j<p;j++)
				tmp += A[j*n + i]*v[k*p + j];
			outVector[k*n+i] = tmp;
		}
	}
}

int MatrixInvSymmetric(double *a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
        }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[(n-1)*n+i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    free(b);
    return(2);
}

double incbeta(double x, double a, double b) {
	if (x < 0.0 || x > 1.0) return 1.0/0.0;

	/*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
	if (x > (a+1.0)/(a+b+2.0)) {
		return (1.0-incbeta(1.0-x,b,a)); /*Use the fact that beta is symmetrical.*/
	}

	/*Find the first part before the continued fraction.*/
	const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
	const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

	/*Use Lentz's algorithm to evaluate the continued fraction.*/
	double f = 1.0, c = 1.0, d = 0.0;

	int i, m;
	for (i = 0; i <= ITMAX; ++i) {
		m = i/2;

		double numerator;
		if (i == 0) {
			numerator = 1.0; /*First numerator is 1.0.*/
		} else if (i % 2 == 0) {
			numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
		} else {
			numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
		}

		/*Do an iteration of Lentz's algorithm.*/
		d = 1.0 + numerator * d;
		if (fabs(d) < FPMIN) d = FPMIN;
		d = 1.0 / d;

		c = 1.0 + numerator / c;
		if (fabs(c) < FPMIN) c = FPMIN;

		const double cd = c*d;
		f *= cd;

		/*Check for stop.*/
		if (fabs(1.0-cd) < EPS) {
			return front * (f-1.0);
		}
	}

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}

//-------------- SLR of GLM regression --------------------------------
int CholeskyDecomL(double *L, int n, double *a){
	// L is lowertriangle matrix satisfying a = L*L'
	int i,j,k;
	double sum, *p;
	p = (double*)malloc(sizeof(double)*n);

	for(i=0; i<n; i++)
	{
		for(j=i; j<n;j++)
		{
			for(sum=a[j*n+i], k=i-1; k>=0; k--)
			{
				sum-=a[k*n+i]*a[k*n+j];
			}
			if (i==j)
			{
				if(sum<=0.0)
					return(-1);
				p[i]=sqrt(sum);
			}
			else
				a[i*n+j]=sum/p[i];
		}
	}
	for(i=0;i<n;i++){
		for(j=i+1; j<n;j++){
			L[i*n+j] = a[i*n+j];
		}
		L[i*n+i] = p[i];
		for(j=0; j<i;j++){
			L[i*n+j] = 0.0;
		}
	}
	free(p);
	return(1);
}

int CholeskyDecomU(double *U, int n, double *a){
	// U is uppertriangle matrix satisfying a = U'*U
	int i,j,k;
	double sum, *p;
	p = (double*)malloc(sizeof(double)*n);

	for(i=0; i<n; i++)
	{
		for(j=i; j<n;j++)
		{
			for(sum=a[j*n+i], k=i-1; k>=0; k--)
			{
				sum-=a[k*n+i]*a[k*n+j];
			}
			if (i==j)
			{
				if(sum<=0.0)
					return(-1);
				p[i]=sqrt(sum);
			}
			else
				a[i*n+j]=sum/p[i];
		}
	}
	for(i=0;i<n;i++){
		for(j=i+1; j<n;j++){
			U[j*n+i] = a[i*n+j];
		}
		U[i*n+i] = p[i];
		for(j=0; j<i;j++){
			U[j*n+i] = 0.0;
		}
	}
	free(p);
	return(1);
}

double GLM_SLR_bstr(double *Tns, double *tx, double *x, double *z, double *resids, double *weight, int *indx,
		double *I0, int n, int p, int p2, int p3, int K, int M, double *G){
	int i,j,k,s,t,p11=p+p2,*subg,sumsb=0;
	double tmp, tmp1, *psi1, *psi0, *C11, *psin, Tn0=0.0;

	subg	= (int*)malloc(sizeof(int)*n);
	psi0 	= (double*)malloc(sizeof(double)*n*p11);
	psi1 	= (double*)malloc(sizeof(double)*n*p11);
	psin 	= (double*)malloc(sizeof(double)*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			psi1[j*n+i]	= tx[j*n+i]*weight[i];
			psi0[j*n+i]	= tx[j*n+i]*resids[i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}
	for(k=0; k<K; k++){
		sumsb 	= 0;
		for(i=0; i<n; i++){
			subg[i] = indx[k*n+i];
			sumsb += subg[i];
		}

		if (sumsb>0){
			for(i=0; i<n; i++){
				if(subg[i]){
					for(j=0; j<p2; j++){
						psi1[(j+p)*n+i]	= x[j*n+i]*weight[i];
						psi0[(j+p)*n+i]	= x[j*n+i]*resids[i];
					}
				}
				else{
					for(j=0; j<p2; j++){
						psi1[(j+p)*n+i]	= 0.0;
						psi0[(j+p)*n+i]	= 0.0;
					}
				}
			}
			for (s = 0; s < p11; s++){
				for (t = s; t < p11; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi1[s*n+i]*psi1[t*n+i];
					}
					C11[s*p11+t] = tmp;
				}
			}
			for (s = 1; s < p11; s++){
				for (t = 0; t < s; t++){
					C11[s*p11+t] = C11[t*p11+s];
				}
			}
			MatrixInvSymmetric(C11,p11);

			for (s = 0; s < p; s++){
				for (t = 0; t < p; t++){
					C11[s*p11+t] -= I0[s*p+t];
				}
			}

			for (s = 0; s < p11; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi0[s*n+i];
				}
				psin[s] = tmp;
			}
			tmp1 = 0.0;
			for (s = 0; s < p11; s++){
				for (t = 0; t < p11; t++){
					tmp1 += psin[s]*C11[s*p11+t]*psin[t];
				}
			}

			if(tmp1 > Tn0){
				Tn0 = tmp1;
			}


			for(j=0; j< M; j++){
				for (s = 0; s < p11; s++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi1[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp;
				}
				tmp1 = 0.0;
				for (s = 0; s < p11; s++){
					for (t = 0; t < p11; t++){
						tmp1 += psin[s]*C11[s*p11+t]*psin[t];
					}
				}

				if(tmp1 > Tns[j]){
					Tns[j] = tmp1;
				}
			}

		}
	}


	free(subg);
	free(psi0);
	free(psi1);
	free(psin);
	free(C11);
	return Tn0;
}

SEXP _GLM_SLR(SEXP tX, SEXP X, SEXP Z, SEXP RESIDS, SEXP WEIGHT, SEXP INDX, SEXP I0, SEXP G, SEXP DIMs){
	int n, p1, p2, p3, K, M;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];

	SEXP rTn0, rTns, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
  	PROTECT(rTns 		= allocVector(REALSXP, 	M));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

	REAL(rTn0)[0] = GLM_SLR_bstr(REAL(rTns), REAL(tX), REAL(X), REAL(Z), REAL(RESIDS), REAL(WEIGHT),
							INTEGER(INDX), REAL(I0), n, p1, p2, p3, K, M, REAL(G));

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("Tnstar"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTns);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- SLR_C of GLM regression --------------------------------
double EstLinearSLR0(double *beta, const double *x, const double *y, double *residual, double *hess, int n, int p){
	int i,j,k;
	double tmp, *xy, loglokl=0.0;
	xy 		= (double*)malloc(sizeof(double)*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		tmp = y[i] - tmp;
		residual[i] = tmp;
		loglokl += tmp*tmp;
	}
	loglokl /= n-p;
	loglokl = log(loglokl) + LOG2PI - n + p;

	free(xy);
	return loglokl;
}

double EstLogisticSLR0(double *beta0, double *x, double *y, double *residual, double *weight, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double tmp, bnorm, yk, expx, wk, loglokl=0.0;
	double *beta, *qy, *dpsi;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
		AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		loglokl += y[i]*tmp - log(1+expx);
		wk 	= expx/(1.0+expx);
		yk 	= y[i] - wk;
		tmp = wk*(1.0-wk);
		residual[i] = yk;
		weight[i] 	= tmp;
		for(j=0;j<p;j++){
			dpsi[j*n+i] = x[j*n+i]*tmp;
		}
	}

	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*dpsi[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}


	MatrixInvSymmetric(hess,p);

	free(beta);
	free(qy);
	free(dpsi);
	return 2.0*loglokl;
}

double EstPoissSLR0(double *beta0, double *x, double *y, double *residual, double *weight, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta, *qy, *dpsi, loglokl=0.0;
	double tmp, bnorm, yk, wk;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			wk 	= exp(tmp);
			yk 	= wk - y[i];
			tmp = wk;
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		wk 	= exp(tmp);
		loglokl += y[i]*tmp - wk;// loglokl += y[i]*tmp - log(1+expx) - log(y[i]!)
		yk 	= y[i] - wk;
		residual[i] = yk;
		weight[i] 	= wk;
		for(j=0;j<p;j++){
			dpsi[j*n+i] = x[j*n+i]*wk;
		}
	}

	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*dpsi[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}
	MatrixInvSymmetric(hess,p);

	free(beta);
	free(qy);
	free(dpsi);
	return 2.0*loglokl;
}

double EstLinearSLR(double *beta, const double *x, const double *y, double *hess, int n, int p){
	int i,j,k;
	double tmp, *xy, loglokl=0.0;
	xy 		= (double*)malloc(sizeof(double)*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		tmp = y[i] - tmp;
		loglokl += tmp*tmp;
	}
	loglokl /= n-p;
	loglokl = log(loglokl) + LOG2PI - n + p;

	free(xy);
	return loglokl;
}

double EstLogisticSLR(double *beta0, double *x, double *y, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double tmp, bnorm, yk, expx, wk, loglokl=0.0;
	double *beta, *qy, *dpsi;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
		AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		loglokl += y[i]*tmp - log(1+expx);
	}

	free(beta);
	free(qy);
	free(dpsi);
	return 2.0*loglokl;
}

double EstPoissSLR(double *beta0, double *x, double *y, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta, *qy, *dpsi, loglokl=0.0;
	double tmp, bnorm, yk, wk;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			wk 	= exp(tmp);
			yk 	= wk - y[i];
			tmp = wk;
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		wk 	= exp(tmp);
		loglokl += y[i]*tmp - wk;// loglokl += y[i]*tmp - log(1+expx) - log(y[i]!)
	}

	free(beta);
	free(qy);
	free(dpsi);
	return 2.0*loglokl;
}

void GLM_SLR_C(double *pvals, double *y, double *tx, double *x, double *z, double *theta, double* G, double *thetahat,
		int n, int p, int p2, int p3, int K, int M, int type, int maxstep, double tol){
	int i,j,k,s,t,p11=p+p2,*subg,sumsb=0,count=0,maxk=0;
	double tmp1, loglokl0, loglokl1, Tn0;
	double *psi0, *psi1, *psi2, *psin, *C11, *alpha0, *weight, *Tns, *I0, *In, *resids;

	subg	= (int*)malloc(sizeof(int)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p11);
	weight 	= (double*)malloc(sizeof(double)*n);
	resids  = (double*)malloc(sizeof(double)*n);
	Tns		= (double*)malloc(sizeof(double)*M);
	psi0 	= (double*)malloc(sizeof(double)*n*p11);
	psi1 	= (double*)malloc(sizeof(double)*n*p11);
	psi2 	= (double*)malloc(sizeof(double)*n*p11);
	psin 	= (double*)malloc(sizeof(double)*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	I0		= (double*)malloc(sizeof(double)*p*p);
	In 		= (double*)malloc(sizeof(double)*p11*p11);

	for(j=0;j<p11;j++)	alpha0[j] 	= 0.0;
	if(type==1){// Gaussian
		loglokl0 = EstLinearSLR0(alpha0, tx, y, resids, I0, n, p);
		for(i=0; i<n; i++){
			weight[i] = 1.0;
		}
	}
	else if(type==2){//binomial
		loglokl0 = EstLogisticSLR0(alpha0, tx, y, resids, weight, I0, n, p, maxstep, tol);
		for(i=0; i<n; i++){
			weight[i] = sqrt(weight[i]);
		}
	}
	else{//Poisson
		loglokl0 = EstPoissSLR0(alpha0, tx, y, resids, weight, I0, n, p, maxstep, tol);
		for(i=0; i<n; i++){
			weight[i] = sqrt(weight[i]);
		}
	}

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			psi0[j*n+i]	= tx[j*n+i]*resids[i];
			psi1[j*n+i]	= tx[j*n+i]*weight[i];
			psi2[j*n+i]	= tx[j*n+i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}
	loglokl1 = -1000000.0;
	Tn0 = -1000000.0;
	for(k=0; k<K; k++){
		sumsb 	= 0;
		if(p3==1){
			for(i=0; i<n; i++){
				subg[i] = IDEX(theta[k], z[i]);
				sumsb += subg[i];
			}
		}
		else{
			for(i=0; i<n; i++){
				tmp1 = 0.0;
				for (s = 0; s < p3; s++){
					tmp1 += z[s*n+i]*theta[k*p3 + s];
				}
				subg[i] = IDEX(0.0, tmp1);
				sumsb += subg[i];
			}
		}

		if (sumsb>0){
			for(i=0; i<n; i++){
				if(subg[i]){
					for(j=0; j<p2; j++){
						psi2[(j+p)*n+i]	= x[j*n+i];
					}
				}
				else{
					for(j=0; j<p2; j++){
						psi2[(j+p)*n+i]	= 0.0;
					}
				}
			}

			if(type==1){// Gaussian
				tmp1 = EstLinearSLR(alpha0, psi2, y, C11, n, p11);
			}
			else if(type==2){//binomial
				tmp1 = EstLogisticSLR(alpha0, psi2, y, C11, n, p11, maxstep, tol);
			}
			else{//Poisson
				tmp1 = EstPoissSLR(alpha0, psi2, y, C11, n, p11, maxstep, tol);
			}
			if(tmp1 > loglokl1){
				loglokl1 = tmp1;
			}

			for(i=0; i<n; i++){
				if(subg[i]){
					for(j=0; j<p2; j++){
						psi0[(j+p)*n+i]	= x[j*n+i]*resids[i];
						psi1[(j+p)*n+i]	= x[j*n+i]*weight[i];
					}
				}
				else{
					for(j=0; j<p2; j++){
						psi0[(j+p)*n+i]	= 0.0;
						psi1[(j+p)*n+i]	= 0.0;
					}
				}
			}

			for (s = 0; s < p11; s++){
				for (t = s; t < p11; t++){
					tmp1 = 0.0;
					for(i=0; i<n; i++){
						tmp1 += psi1[s*n+i]*psi1[t*n+i];
					}
					In[s*p11+t] = tmp1;
				}
			}
			for (s = 1; s < p11; s++){
				for (t = 0; t < s; t++){
					In[s*p11+t] = In[t*p11+s];
				}
			}

			MatrixInvSymmetric(In, p11);

			for (s = 0; s < p; s++){
				for (t = 0; t < p; t++){
					In[s*p11+t] -= I0[s*p+t];
				}
			}

			for (s = 0; s < p11; s++){
				tmp1 = 0.0;
				for(i=0; i<n; i++){
					tmp1 += psi0[s*n+i];
				}
				psin[s] = tmp1;
			}
			tmp1 = 0.0;
			for (s = 0; s < p11; s++){
				for (t = 0; t < p11; t++){
					tmp1 += psin[s]*In[s*p11+t]*psin[t];
				}
			}

			if(tmp1 > Tn0){
				Tn0 = tmp1;
				maxk = k;
			}

			for(j=0; j< M; j++){
				for (s = 0; s < p11; s++){
					tmp1 = 0.0;
					for(i=0; i<n; i++){
						tmp1 += psi1[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp1;
				}
				tmp1 = 0.0;
				for (s = 0; s < p11; s++){
					for (t = 0; t < p11; t++){
						tmp1 += psin[s]*In[s*p11+t]*psin[t];
					}
				}

				if(tmp1 > Tns[j]){
					Tns[j] = tmp1;
				}
			}

		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals[0] = 1.0*count/M;

	count = 0;
	loglokl1 -= loglokl0;
	for(j=0; j< M; j++){
		if(loglokl1<=Tns[j]){
			count++;
		}
	}
	pvals[1] = 1.0*count/M;
	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}

	free(subg);
	free(alpha0);
	free(weight);
	free(resids);
	free(Tns);
	free(psi0);
	free(psi1);
	free(psi2);
	free(psin);
	free(C11);
	free(I0);
	free(In);
}

void GLM_SLR_Approx_C(double *pvals, double *y, double *tx, double *x, double *z, double *theta, double* G, double *thetahat,
		int n, int p, int p2, int p3, int K, int M, int type, int maxstep, double tol){
	int i,j,k,s,t,p11=p+p2,*subg,sumsb=0,count=0,maxk=0;
	double tmp1, Tn0;
	double *psi0, *psi1, *psi2, *psin, *C11, *alpha0, *weight, *Tns, *I0, *In, *resids;

	subg	= (int*)malloc(sizeof(int)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p11);
	weight 	= (double*)malloc(sizeof(double)*n);
	resids  = (double*)malloc(sizeof(double)*n);
	Tns		= (double*)malloc(sizeof(double)*M);
	psi0 	= (double*)malloc(sizeof(double)*n*p11);
	psi1 	= (double*)malloc(sizeof(double)*n*p11);
	psi2 	= (double*)malloc(sizeof(double)*n*p11);
	psin 	= (double*)malloc(sizeof(double)*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	I0		= (double*)malloc(sizeof(double)*p*p);
	In 		= (double*)malloc(sizeof(double)*p11*p11);

	for(j=0;j<p11;j++)	alpha0[j] 	= 0.0;
	if(type==1){// Gaussian
		EstLinearSLR0(alpha0, tx, y, resids, I0, n, p);
		for(i=0; i<n; i++){
			weight[i] = 1.0;
		}
	}
	else if(type==2){//binomial
		EstLogisticSLR0(alpha0, tx, y, resids, weight, I0, n, p, maxstep, tol);
		for(i=0; i<n; i++){
			weight[i] = sqrt(weight[i]);
		}
	}
	else if(type==3){//Poisson
		EstPoissSLR0(alpha0, tx, y, resids, weight, I0, n, p, maxstep, tol);
		for(i=0; i<n; i++){
			weight[i] = sqrt(weight[i]);
		}
	}

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			psi0[j*n+i]	= tx[j*n+i]*resids[i];
			psi1[j*n+i]	= tx[j*n+i]*weight[i];
			psi2[j*n+i]	= tx[j*n+i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}
	Tn0 = -1000000.0;
	for(k=0; k<K; k++){
		sumsb 	= 0;
		if(p3==1){
			for(i=0; i<n; i++){
				subg[i] = IDEX(theta[k], z[i]);
				sumsb += subg[i];
			}
		}
		else{
			for(i=0; i<n; i++){
				tmp1 = 0.0;
				for (s = 0; s < p3; s++){
					tmp1 += z[s*n+i]*theta[k*p3 + s];
				}
				subg[i] = IDEX(0.0, tmp1);
				sumsb += subg[i];
			}
		}

		if (sumsb>0){
			for(i=0; i<n; i++){
				if(subg[i]){
					for(j=0; j<p2; j++){
						psi2[(j+p)*n+i]	= x[j*n+i];
					}
				}
				else{
					for(j=0; j<p2; j++){
						psi2[(j+p)*n+i]	= 0.0;
					}
				}
			}

			for(i=0; i<n; i++){
				if(subg[i]){
					for(j=0; j<p2; j++){
						psi0[(j+p)*n+i]	= x[j*n+i]*resids[i];
						psi1[(j+p)*n+i]	= x[j*n+i]*weight[i];
					}
				}
				else{
					for(j=0; j<p2; j++){
						psi0[(j+p)*n+i]	= 0.0;
						psi1[(j+p)*n+i]	= 0.0;
					}
				}
			}

			for (s = 0; s < p11; s++){
				for (t = s; t < p11; t++){
					tmp1 = 0.0;
					for(i=0; i<n; i++){
						tmp1 += psi1[s*n+i]*psi1[t*n+i];
					}
					In[s*p11+t] = tmp1;
				}
			}
			for (s = 1; s < p11; s++){
				for (t = 0; t < s; t++){
					In[s*p11+t] = In[t*p11+s];
				}
			}


			MatrixInvSymmetric(In, p11);

			for (s = 0; s < p; s++){
				for (t = 0; t < p; t++){
					In[s*p11+t] -= I0[s*p+t];
				}
			}

			for (s = 0; s < p11; s++){
				tmp1 = 0.0;
				for(i=0; i<n; i++){
					tmp1 += psi0[s*n+i];
				}
				psin[s] = tmp1;
			}
			tmp1 = 0.0;
			for (s = 0; s < p11; s++){
				for (t = 0; t < p11; t++){
					tmp1 += psin[s]*In[s*p11+t]*psin[t];
				}
			}

			if(tmp1 > Tn0){
				Tn0 = tmp1;
				maxk = k;
			}

			for(j=0; j< M; j++){
				for (s = 0; s < p11; s++){
					tmp1 = 0.0;
					for(i=0; i<n; i++){
						tmp1 += psi1[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp1;
				}
				tmp1 = 0.0;
				for (s = 0; s < p11; s++){
					for (t = 0; t < p11; t++){
						tmp1 += psin[s]*In[s*p11+t]*psin[t];
					}
				}

				if(tmp1 > Tns[j]){
					Tns[j] = tmp1;
				}
			}

		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals[0] = 1.0*count/M;

	pvals[1] = pvals[0];
	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}

	free(subg);
	free(alpha0);
	free(weight);
	free(resids);
	free(Tns);
	free(psi0);
	free(psi1);
	free(psi2);
	free(psin);
	free(C11);
	free(I0);
	free(In);
}

SEXP _GLM_SLR_C(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP THETA, SEXP G, SEXP DIMs, SEXP TOL){
	int n, p1, p2, p3, K, M, type, maxstep, isApprox;
	double tol;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	type 	= INTEGER(DIMs)[6];
	maxstep = INTEGER(DIMs)[7];
	isApprox= INTEGER(DIMs)[8];
	tol 	= REAL(TOL)[0];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	2));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

	if(isApprox){
		GLM_SLR_Approx_C(REAL(rpvals), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(THETA), REAL(G), REAL(rthetahat),
				n, p1, p2, p3, K, M, type, maxstep, tol);
	}
	else{
		GLM_SLR_C(REAL(rpvals), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(THETA), REAL(G), REAL(rthetahat),
				n, p1, p2, p3, K, M, type, maxstep, tol);
	}
	SET_STRING_ELT(list_names, 	0,	mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- WAST of GLM regression --------------------------------
void EstLinearWAST00(double *beta, const double *x, const double *y, double *residual, double *hess, int n, int p){
	int i,j,k;
	double tmp, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(xy);
}

void EstLinearWAST10(const double *x, double *hess, int n, int p){
	int i,j,k;
	double tmp;
	for(j=0; j < p; j++){
		for(k=j; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	for(j=1; j < p; j++){
		for(k=0; k < j; k++){
			hess[j*p+k] = hess[k*p+j];
		}
	}
	MatrixInvSymmetric(hess,p);
}

void EstLinearWAST1(double *beta, const double *x, const double *y, double *residual, double *hess, int n, int p){
	int i,j;
	double tmp, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(xy);
}

void EstLogisticWAST1(double *beta0, double *x, double *y, double *residual, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double tmp, bnorm, yk, expx, wk;
	double *beta, *qy, *dpsi;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
		AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		wk 	= expx/(1.0+expx);
		residual[i] = y[i] - wk;
	}

	free(beta);
	free(qy);
	free(dpsi);
}

void EstPoissWAST1(double *beta0, double *x, double *y, double *residual, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta, *qy, *dpsi;
	double tmp, bnorm, yk, wk;

	beta 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			wk 	= exp(tmp);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*wk;
			}
		}
		for(j=0; j < p; j++){
			for(k=j; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		for(j=1; j < p; j++){
			for(k=0; k < j; k++){
				hess[j*p+k] = hess[k*p+j];
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - exp(tmp);
	}

	free(beta);
	free(qy);
	free(dpsi);
}

double GLM_single_wast(double *Tns, double *y, double *tx, double *x, double *z, double *resid, int n, int p, int p2, int p3,
		int isBeta, double shape1, double shape2, int maxIter, double tol, int type, int M){
	int i,j,s,g;
	double tmp, omega, Tn, Tn0=0.0, xij;
	double *ty, *C11, *OMEGA, *yb, *alpha0;

	yb 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	C11 	= (double*)malloc(sizeof(double)*p*p);
	ty  	= (double*)malloc(sizeof(double)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);


	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf(SQRT2*z[i]) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			tmp = omega*xij;
			OMEGA[i*n+j]	= tmp;
			Tn0 	+= tmp*resid[i]*resid[j];
		}
	}
	if(type==1){// Gaussian
		EstLinearWAST10(tx, C11, n, p);
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}

		if(type==2){//binomial
			for(j=0;j<p;j++)	alpha0[j] 	= 0.0;
		}
		else if(type==3){//Poisson
			for(j=0;j<p;j++)	alpha0[j] 	= 0.0;
		}

		if(type==2){//binomial
			EstLogisticWAST1(alpha0, tx, yb, resid, C11, n, p, maxIter, tol);
		}
		else if(type==3){//Poisson
			EstPoissWAST1(alpha0, tx, yb, resid, C11, n, p, maxIter, tol);
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(alpha0);
	free(ty);
	free(OMEGA);
	free(C11);
	return Tn0;
}

double GLM_multiple_wast(double *Tns, double *y, double *tx, double *x, double *z, double *resid, int n, int p,
		int p2, int p3, double shape1, double shape2, int maxIter, double tol, int type, int M){
	int i,j,s,g;
	double tmp, omega, Tn, Tn0=0.0, xij, sd, rho;
	double *stdx, *OMEGA, *yb, *alpha0, *C11;

	yb 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);


	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			tmp = omega*xij;
			OMEGA[i*n+j]	= tmp;
			Tn0 	+= tmp*resid[i]*resid[j];
		}
	}

	if(type==1){// Gaussian
		EstLinearWAST10(tx, C11, n, p);
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}

		if(type==2){//binomial
			for(j=0;j<p;j++)	alpha0[j] 	= 0.0;
		}
		else if(type==3){//Poisson
			for(j=0;j<p;j++)	alpha0[j] 	= 0.0;
		}

		if(type==1){// Gaussian
			EstLinearWAST1(alpha0, tx, yb, resid, C11, n, p);
		}
		else if(type==2){//binomial
			EstLogisticWAST1(alpha0, tx, yb, resid, C11, n, p, maxIter, tol);
		}
		else if(type==3){//Poisson
			EstPoissWAST1(alpha0, tx, yb, resid, C11, n, p, maxIter, tol);
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}
		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(alpha0);
	free(stdx);
	free(OMEGA);
	free(C11);
	return 2.0*Tn0/n;
}

SEXP _GLM_WAST(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESIDS, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxIter, isBeta, type, M;
	double shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];
	type 	= INTEGER(DIMs)[7];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn0, rTns, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTns 		= allocVector(REALSXP, 	M));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	if(p3==1){
		REAL(rTn0)[0] = GLM_single_wast(REAL(rTns), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESIDS),
						n, p1, p2, p3, isBeta, shape1, shape2, maxIter, tol, type, M);
	}
	else{
		REAL(rTn0)[0] = GLM_multiple_wast(REAL(rTns), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESIDS),
						n, p1, p2, p3, shape1, shape2, maxIter, tol, type, M);
	}

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("Tns"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTns);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}
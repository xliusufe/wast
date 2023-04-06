#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <float.h>  	// required for DBL_EPSILON
#include <math.h>
#include "_WAST_HEAD.h"

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

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q){
    int i,j,k;
    double temp;
	for (i = 0; i<p; i++){
		for (k = 0; k<q; k++){
			temp = 0.0;
			for (j = 0; j < n; j++)
				temp += A[i + j*p] * B[k + j*q];
			outMatrix[i*q + k] = temp;
		}
	}
}

double exp_approx(double x, int n){
	x = 1.0 + x/256;
	for (int i = 0; i < n; i++){
		x *= x;
	}

	return x;
}

void Standarize(double *y, double *std, int n, int p, double *x, int flag)
{
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) s += x[j*n+i];
		s = s/n;
		for(i=0;i<n;i++) s1  += x[j*n+i]*x[j*n+i];
		s1 = s1/n - s*s;
		if(flag)
			std[j] = sqrt(s1);
		else
			std[j] = sqrt(n*s1/(n-1));
		for(i=0;i<n;i++) y[j*n+i] = (x[j*n+i]-s)/std[j];
	}
}

void Std(double *std, double *x, int n, int p, int flag){
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) s += x[j*n+i];
		s = s/n;
		for(i=0;i<n;i++) s1  += x[j*n+i]*x[j*n+i];
		s1 = s1/n - s*s;
		if(flag)
			std[j] = sqrt(s1);
		else
			std[j] = sqrt(n*s1/(n-1));
	}
}

void sortN(int *ind0, double *x, int n, int dd){
	int i, j, MaxInd, d, *ind;
	double tmp;
	ind = (int*)malloc(sizeof(int)*n);
	for(i=0;i<n;i++) ind[i] = i;

	d = (dd==n?dd-1:dd);
	for(i=0;i<d;i++)
	{
		tmp = x[0]; MaxInd = ind[0];
		for(j=1;j<n-i;j++)
		{
			if(x[j]<tmp)
			{
				x[j-1] = x[j];
				x[j] = tmp;
				ind[j-1] = ind[j];
				ind[j] = MaxInd;
			}
			else
			{
				tmp = x[j];
				MaxInd = ind[j];
			}
		}
	}
	for(j=0;j<dd;j++) ind0[j] = ind[n-j-1];
	free(ind);
}

void SortQ(double *s, int l, int r)
{
    int i, j;
	double x;
    if (l < r)
    {
        i = l;
        j = r;
        x = s[i];
        while (i < j)
        {
            while(i < j && s[j] > x) j--;
			if(i < j) s[i++] = s[j];
            while(i < j && s[i] < x) i++;
			if(i < j) s[j--] = s[i];
        }
        s[i] = x;
        SortQ(s, l, i-1);
        SortQ(s, i+1, r);
    }
}

int LowTriangularInv(double *B, int n, double *A){
	// Input:
	// A is a lower triangular matrix
	//
	// Output:
	// B = inv(A)
	//
	int i,j,k;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){
	// Input:
	// X is a p*n matrix
	//
	// Output:
	// R is a p*p lower triangular matrix
	// E is a p*n matrix satisfying E*t(E) = I_p
	//
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;

	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p
	for(j=1;j<p;j++){
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

void SampleQuantile(double *qr, int m, double *z, int n, double *q)
{
	double *zs=(double*)malloc(sizeof(double)*n);
	int i, ind;
	for(i=0;i<n;i++) zs[i] = z[i];
	SortQ(zs, 0, n-1);
	for(i=0;i<m;i++)
	{
		ind = floor(q[i]*n);
		if (ind!=n*q[i])
			qr[i] = zs[ind];
		else
			qr[i] = (zs[ind-1] + zs[ind])/2;
	}
	free(zs);
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

	// CholeskyDecomU(H, p, I0);
	// for (s = 0; s < p; s++){
	// 	for (t = 0; t < p; t++){
	// 		tmp = 0.0;
	// 		for(j=0; j< t+1; j++){
	// 			tmp += H[s*p+j]*H[t*p+j];
	// 		}
	// 		printf("%f  ", tmp);
	// 	}
	// 	printf("\n");
	// }
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

void Kernelh(double *x, double *kern, int n, double h0, int type){
	// type = 1: Gaussian kernel
	// type = 2; Epanecknikov kernel
	// type = 3: uniform kernel
	int i;
	double tmp, h;
	Std( &h, x, n, 1, 0);
	h /= h0;
	for(i=0; i<n; i++){
		tmp = x[i]*h;
		tmp *= tmp;
		if(type==1){
			kern[i] = exp(-0.5*tmp)*MPI2*h;
		}
		else if(type==2){
			kern[i] = 0.75*(1-tmp)*(tmp<1?1:0)*h;
		}
		else{
			kern[i] = (tmp>0.5?0:1)*(tmp<-0.5?0:1)*h;
		}
	}
}

void EstQR2(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, q2 = 2*qq-1, yk, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	SampleQuantile(&yk, 1, y, n, &qq);
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(i=0;i<n;i++)	residual[i]	= y[i] - yk;

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

	if(p<2){
		if(hess[0]<MEPS){
			step = maxstep;
		}
		else{
			hess[0] = 1.0/hess[0];
		}
	}
	else{
		flag = MatrixInvSymmetric(hess,p);
		if(flag<0){
			step = maxstep;
		}
	}

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			wk = EPS + fabs(residual[i]);
			yk = y[i] + wk*q2;
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
			}
		}

		if(p<2){
			beta[0] = qy[0]*hess[0];
		}
		else{
			AbyB(beta, hess, qy, p, p, 1);
		}

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
				residual[i] = y[i] - tmp;
			}
		}
	}

	free(beta0);
	free(dpsi);
	free(hess);
	free(qy);
}

void EstLogisticR(double *beta, double *residual, double *x, double *y, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;

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
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*wk*(1.0-wk);
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
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
		yk 	= y[i] - wk;
		residual[i] = yk;
	}

	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
}

void EstLinearR(double *beta, double *residual, const double *x, const double *y, int n, int p){
	int i,j,k;
	double tmp, *hess, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0; j < p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(hess);
	free(xy);
}

double EstQuantileR1(double *y, double *residual, double tau, const int n){
	double beta=0.0;
	SampleQuantile(&beta, 1, y, n, &tau);

	for(int i=0;i<n;i++){
		residual[i] = y[i] - beta;
	}

	return beta;
}

void EstQuantileR(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
	double *beta0, *hess, *qy;
	double tmp, bnorm, q2 = 2*qq-1, yk, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);

	SampleQuantile(&yk, 1, y, n, &qq);
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(i=0;i<n;i++)	residual[i]	= y[i] - yk;

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

	if(p<2){
		if(hess[0]<MEPS){
			step = maxstep;
		}
		else{
			hess[0] = 1.0/hess[0];
		}
	}
	else{
		flag = MatrixInvSymmetric(hess,p);
		if(flag<0){
			step = maxstep;
		}
	}

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			wk = EPS + fabs(residual[i]);
			yk = y[i] + wk*q2;
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
			}
		}

		if(p<2){
			beta[0] = qy[0]*hess[0];
		}
		else{
			AbyB(beta, hess, qy, p, p, 1);
		}

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
				residual[i] = y[i] - tmp;
			}
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
		residual[i] = y[i] - tmp;
	}

	free(beta0);
	free(hess);
	free(qy);
}

void Est_probitR(int *y, double *x, double *alpha, double *residual, int n, int p, int maxIter, double tol){
	int i,j, k,step=0;
	double tmp, tmp1, phix1, phix2, bnorm, bnorm0 = 1.0;
	double *mu, *alpha0, *hess, *dpsi, *psi10;

	mu 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	psi10 	= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alpha0[j] = 1.0;
	}
	while(step<maxIter){
		step++;
		for(j=0; j < p; j++){
			psi10[j] = 0.0;
		}
		AbyB(mu, x, alpha0, n, p, 1);
		for(i=0; i<n; i++){
			tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
			tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
			phix1 	= tmp/tmp1;
			phix2 	= tmp/(1-tmp1);
			tmp 	= y[i]?phix1:(-phix2);
			tmp1 	= y[i]?(-phix1*(mu[i]+phix1)):(phix2*(mu[i] - phix2));
			for(j=0; j < p; j++){
				dpsi[j*n+i] = tmp1*x[j*n+i];
				psi10[j] 	+= tmp*x[j*n+i];
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		MatrixInvSymmetric(hess,p);
    	AbyB(alpha, hess, psi10, p, p, 1);

		bnorm = 0.0;
		for(j=0; j < p; j++){
			alpha[j]	= alpha0[j] - alpha[j];
			bnorm 		+= alpha[j]*alpha[j];
		}

		if(sqrt(bnorm/bnorm0)<tol){
			break;
		}
		else{
			bnorm0 = bnorm;
			for(j=0; j < p; j++){
				alpha0[j] = alpha[j];
			}
		}
	}

	AbyB(mu, x, alpha, n, p, 1);
	for(i=0; i<n; i++){
		tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
		tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
		phix1 	= tmp/tmp1;
		phix2 	= - tmp/(1-tmp1);
		tmp 	= y[i]?phix1:phix2;
		residual[i] = tmp;
	}

	free(mu);
	free(alpha0);
	free(psi10);
	free(dpsi);
	free(hess);
}

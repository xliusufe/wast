#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "_WAST_HEAD.h"

//-------------- SST of GLM regression --------------------------------
double SST_GLM_bstr(double *y, double *tx1, double *x, double *z, double *resid1, double *weight, double* thetahat,
		int n, int p11, int p2, int p3, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, tmp2, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *score1, *psi, *psin;
	double *C11, *psi_h, *K1, *Vh, *Tns;

	subg	= (int*)malloc(sizeof(int)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);


	for(j=0; j < p11; j++){
		for(k=0; k < p11; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += tx1[j*n+i]*weight[i]*tx1[k*n+i];
			}
			C11[j*p11+k] = tmp/n;
		}
	}
	MatrixInvSymmetric(C11,p11);

	for(i=0; i<n; i++){
		for(j=0; j<p11; j++){
			tmp = 0.0;
			for (s = 0; s < p11; s++){
				tmp += C11[j*p11+s]*tx1[s*n+i];
			}
			score1[j*n+i] = tmp*resid1[i];
		}

		for(j=0; j<p2; j++){
			psi[j*n+i] 	= x[j*n+i]*resid1[i];
		}
	}
	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

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
				tmp = 0.0;
				for (j = 0; j < p3; j++){
					tmp += z[j*n+i]*theta[k*p3 + j];
				}
				subg[i] = IDEX(0.0, tmp);
				sumsb += subg[i];
			}
		}

		if (sumsb>0){
			for (s = 0; s < p2*p11; s++){
				K1[s] = 0.0;
			}
			for(i=0; i<n; i++){
				if(subg[i]){
					for (s = 0; s < p2; s++){
						for (t = 0; t < p11; t++){
							K1[s*p11+t] += x[s*n+i]*tx1[t*n+i]*weight[i];
						}
					}
				}
			}
			for (s = 0; s < p2*p11; s++){
				K1[s] /= n;
			}

			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp2 = (subg[i])?psi[s*n+i]:0.0;
					tmp += tmp2;
					tmp1 = 0.0;
					for (t = 0; t < p11; t++){
						tmp1 += K1[s*p11+t]*score1[t*n+i];
					}
					psi_h[s*n+i] = tmp2 - tmp1;

				}
				psin[s] = tmp;
			}

			for (s = 0; s < p2; s++){
				for (t = s; t < p2; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi_h[s*n+i]*psi_h[t*n+i];
					}
					Vh[s*p2+t] = tmp/n;
				}
			}
			for (s = 1; s < p2; s++){
				for (t = 0; t < s; t++){
					Vh[s*p2+t] = Vh[t*p2+s];
				}
			}

			if(p2 ==1){
				Vh[0] = 1.0/Vh[0];
			}
			else{
				MatrixInvSymmetric(Vh, p2);
			}

			Tn = 0.0;
			for (s = 0; s < p2; s++){
				for (t = 0; t < p2; t++){
					Tn += psin[s]*Vh[s*p2+t]*psin[t];
				}
			}

			if(Tn>Tn0){
				Tn0 = Tn;
				maxk = k;
			}
			for(j=0; j< M; j++){
				for (s = 0; s < p2; s++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi_h[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp;
				}

				Tn_star = 0.0;
				for (s = 0; s < p2; s++){
					for (t = 0; t < p2; t++){
						Tn_star += psin[s]*Vh[s*p2+t]*psin[t];
					}
				}
				if(Tn_star>Tns[j]){
					Tns[j] = Tn_star;
				}
			}
		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}

	free(subg);
	free(psi);
	free(score1);
	free(C11);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(Tns);
	return pvals;
}

SEXP _SST_GLM_bstr(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESID, SEXP WEIGHT, SEXP THETA, SEXP G, SEXP DIMs){
	int n, p1, p2, p3, K, M;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	REAL(rpvals)[0] = SST_GLM_bstr(REAL(Y), REAL(tX),
						REAL(X), REAL(Z), REAL(RESID), REAL(WEIGHT), REAL(rthetahat),
						n, p1, p2, p3, K, M, REAL(G), REAL(THETA));



	SET_STRING_ELT(list_names, 	0,	mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- SLR of GLM regression --------------------------------

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
	if(type==1){
		tmp = 0.0;
		for (i = 0; i < n; i++){
			tmp += resid[i]*resid[i];
		}
		Tn0 *= 1.0*(n-p);
		Tn0 /= tmp*tmp;
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

		if(type==1){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += resid[i]*resid[i];
			}
			Tns[g] *= 1.0*(n-p);
			Tns[g] /= tmp*tmp;
		}

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
	if(type==1){
		tmp = 0.0;
		for (i = 0; i < n; i++){
			tmp += resid[i]*resid[i];
		}
		Tn0 *= 1.0*(n-p);
		Tn0 /= tmp*tmp;
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

		if(type==1){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += resid[i]*resid[i];
			}
			Tns[g] *= 1.0*(n-p);
			Tns[g] /= tmp*tmp;
		}
	}

	free(yb);
	free(alpha0);
	free(stdx);
	free(OMEGA);
	free(C11);
	return 2.0*Tn0/n;
}

double GLM_multiple_wast_approx(double *Tns, double *y, double *tx, double *x, double *z, double *resid, double *zk, double *mu0, int n, int p,
		int p2, int p3, double shape1, double shape2, int maxIter, double tol, int type, int typeGam, int M, int N0){
	int i,j,s,g,count;
	double tmp, tmp1, tmp2, tmp3, omega=0.0, Tn, Tn0=0.0, xij, sd, rho;
	double *stdx, *OMEGA, *yb, *alpha0, *C11, *zmu;

	yb 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	zmu 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);


	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		tmp = 0.0;
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
			tmp += stdx[j*n+i]*mu0[j];
		}
		zmu[i] = tmp;
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
				tmp1 = zmu[j];
				count = 0;
				for (s = 0; s < N0; s++) {
					if(zk[s] < tmp1)
						count++;
				}
				omega = 1.0*count/N0;
			}
			else{
				tmp		= 0.0;
				tmp1 	= zmu[j];
				tmp2 	= rho/sqrt(1-rho*rho);
				tmp3 	= tmp1/sqrt(1-rho*rho);
				for (s = 0; s < N0; s++) {
					if(zk[s] < tmp1){
						tmp += erf(SQRT2*(tmp3 - zk[s]*tmp2)) + 1.0;
					}
				}
				omega = 0.5*tmp/N0;
			}
			tmp = omega*xij;
			OMEGA[i*n+j]	= tmp;
			Tn0 	+= tmp*resid[i]*resid[j];
		}
	}

	if(type==1){
		tmp = 0.0;
		for (i = 0; i < n; i++){
			tmp += resid[i]*resid[i];
		}
		Tn0 *= 1.0*(n-p);
		Tn0 /= tmp*tmp;
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
		else{//Poisson
			EstPoissWAST1(alpha0, tx, yb, resid, C11, n, p, maxIter, tol);
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}
		Tns[g] = 2.0*Tn/n;

		if(type==1){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += resid[i]*resid[i];
			}
			Tns[g] *= 1.0*(n-p);
			// Tns[g] /= tmp;
			Tns[g] /= tmp*tmp;
		}

	}

	free(yb);
	free(alpha0);
	free(zmu);
	free(stdx);
	free(OMEGA);
	free(C11);
	return 2.0*Tn0/n;
}

SEXP _GLM_WAST_APPROX(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESIDS, SEXP Z_K, SEXP MU0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxIter, isBeta, type, typeGam, M, N0;
	double shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];
	type 	= INTEGER(DIMs)[7];
	typeGam = INTEGER(DIMs)[8];
	N0 		= INTEGER(DIMs)[9];

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
		REAL(rTn0)[0] = GLM_multiple_wast_approx(REAL(rTns), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESIDS), REAL(Z_K), REAL(MU0),
						n, p1, p2, p3, shape1, shape2, maxIter, tol, type, typeGam, M, N0);
	}

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("Tns"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTns);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
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

//-------------- ESTIMATION of GLM regression --------------------------------
void EstLogisticCP(double *beta0, double *x, double *y, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
	double tmp, bnorm, yk, wk;
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
			wk 	= 1.0/(1.0+exp(-tmp));
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

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = qy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, qy, p, p, 1);
			}
		}

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

	free(beta);
	free(qy);
	free(dpsi);
}

void EstPoissCP(double *beta0, double *x, double *y, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
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

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = qy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, qy, p, p, 1);
			}
		}

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

	free(beta);
	free(qy);
	free(dpsi);
}

void _EstLinearCP(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *xy, *zy;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}
		

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}


	while(step < maxstep){
		step++;
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			yk 	= y[i] - tmp;
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*wk*yk;
				wz[j*n+i] = z[n+j*n+i]*wk;
			}
			ologelr += yk*yk;
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += wz[j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			nlogelr += tmp*tmp;
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
}

void _EstLogisticCP(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, yk, expx, wk, ologelr=10000*n, nlogelr;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *zy;

	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	while (step < maxstep){
		step++;

		EstLogisticCP(beta0, w, y, hess, n, p, maxstep, eps);

		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk 	= 1.0/(1.0+exp(-tmp));

			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk 	= y[i] - wk;
			tmp = tmp*dsh[i];
			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*tmp*yk;
				wz[j*n+i] = z[n+j*n+i]*tmp*tmp*wk*(1.0-wk);
			}
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += z[n+j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			nlogelr += y[i]*tmp - log(1+expx);
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(zy);
}

void _EstPoissCP(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, yk, wk, ologelr=10000*n, nlogelr;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *zy;

	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	while (step < maxstep){
		step++;

		EstPoissCP(beta0, w, y, hess, n, p, maxstep, eps);

		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk = exp(tmp);
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk 	= y[i] - wk;
			tmp = tmp*dsh[i];
			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*tmp*yk;
				wz[j*n+i] = z[n+j*n+i]*tmp*tmp*wk;
			}
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += z[n+j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk = exp(tmp);
			nlogelr += y[i]*tmp - wk;// loglokl += y[i]*tmp - log(1+expx) - log(y[i]!)
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(zy);
}

SEXP _EST_LINEAR(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	_EstLinearCP(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_LOGISTICR(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	_EstLogisticCP(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_POISSON(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	_EstPoissCP(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- ESTIMATION of GLM regression by Bootstrap method --------------------------------
void EstLogisticCP_Weight(double *beta0, double *x, double *y, double *g, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
	double tmp, bnorm, yk, wk;
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
			wk 	= 1.0/(1.0+exp(-tmp));
			yk 	= wk - y[i];
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk*g[i];
				dpsi[j*n+i] = x[j*n+i]*tmp*g[i];
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

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = qy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, qy, p, p, 1);
			}
		}

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

	free(beta);
	free(qy);
	free(dpsi);
}

void EstPoissCP_Weight(double *beta0, double *x, double *y, double *g, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
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
			yk 	*= g[i];
			wk 	*= g[i];
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

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = qy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, qy, p, p, 1);
			}
		}

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

	free(beta);
	free(qy);
	free(dpsi);
}

void _EstLinearCP_Boot(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *xy, *zy;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp*g[i];
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}


	while(step < maxstep){
		step++;
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = zy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			yk = y[i];
			for(j=0;j<p1;j++){
				yk -= tx[j*n+i]*beta0[j];
			}

			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*wk*yk*g[i];
				wz[j*n+i] = z[n+j*n+i]*wk;
			}
			ologelr += yk*yk*g[i];
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += wz[j*n+i]*wz[k*n+i]*g[i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk = z[i];
			for(j=0;j<p3;j++){
				yk	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				yk 		= 1.0/(1.0+exp(-yk*h));
				sh[i] 	= yk;
				dsh[i] 	= h*yk*(1-yk);
			}
			else if(smooth==2){
				tmp1 	= yk*h;
				yk 		= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= yk;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= yk*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				yk 		= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= yk;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*yk*g[i];
				w[(j+p1)*n+i] = x[j*n+i]*yk;
			}
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			nlogelr += tmp*tmp*g[i];
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
}

void _EstLogisticCP_Boot(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, yk, expx, wk, ologelr=10000*n, nlogelr;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *zy;

	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	while (step < maxstep){
		step++;

		EstLogisticCP_Weight(beta0, w, y, g, hess, n, p, maxstep, eps);

		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk 	= 1.0/(1.0+exp(-tmp));

			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk 	= y[i] - wk;
			yk 	*= g[i];
			tmp = tmp*dsh[i];
			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*tmp*yk;
				wz[j*n+i] = z[n+j*n+i]*tmp*tmp*wk*(1.0-wk)*g[i];
			}
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += z[n+j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			nlogelr += (y[i]*tmp - log(1+expx))*g[i];
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(zy);
}

void _EstPoissCP_Boot(double *beta0, double *theta0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, yk, wk, ologelr=10000*n, nlogelr;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *zy;

	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	while (step < maxstep){
		step++;

		EstPoissCP_Weight(beta0, w, y, g, hess, n, p, maxstep, eps);

		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk = exp(tmp);
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk 	= y[i] - wk;
			yk 	*= g[i];
			wk 	*= g[i];
			tmp = tmp*dsh[i];
			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*tmp*yk;
				wz[j*n+i] = z[n+j*n+i]*tmp*tmp*wk;
			}
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += z[n+j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) + 0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) + 0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += w[j*n+i]*beta0[j];
			}
			wk = exp(tmp);
			nlogelr += (y[i]*tmp - wk)*g[i];// loglokl += y[i]*tmp - log(1+expx) - log(y[i]!)
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(zy);
}

SEXP _EST_LINEAR_Boot(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	_EstLinearCP_Boot(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_LOGISTICR_Boot(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	_EstLogisticCP_Boot(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_POISSON_Boot(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	_EstPoissCP_Boot(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}


//-------------- ESTIMATION of GLM regression with multiple change-planes--------------------------------
void _EstLinearCP_Multi0(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *xy, *zy, *ay, *resids, *resids1, *xb;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	wa 		= (double*)malloc(sizeof(double)*n*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	// for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	// for(j=0;j<p3;j++)	theta0[j] 	= 1.0;
	// for(j=0;j<ng;j++)	ha0[j] 		= (1.0*j)/ng;
	// ha0[0] = -0.5; ha0[1] = 0.5;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	// printf("p = %d, p3 = %d\n", p, p3);

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			

			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// printArrayDouble(xy, p, p);
	// printArrayDouble1(w, n, 5, p);
	while(step < maxstep){
		step++;

		// calculate alpha and beta
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		// printArrayDouble(beta0, p, p);

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k*n+i] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			ologelr += yk*yk;
			for(k=0; k < ng; k++){
				ay[k] += wa[k*n+i]*yk;
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j*n+i]*wa[k*n+i];
				}
			}			
		}

		// printArrayDouble1(hessa, ng, ng, ng);
		// printArrayDouble(ay, ng, ng);

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}
		// printArrayDouble(ha, ng, ng);

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			// printf("nlogelr = %f,   ologelr = %f\n", nlogelr, ologelr);
			// printArrayDouble(resids1, 10, 10);
			// if(step<2){
			// 	printf("nlogelr = %f,   ologelr = %f\n", nlogelr, ologelr);
			// }
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}


		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*wk*yk;
				wz[j*n+i] = z[n+j*n+i]*wk;
			}
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += wz[j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}
		// printArrayDouble1(wz, n, 5, p3);
		// printArrayDouble(zy, p3, p3);
		// printArrayDouble1(hessz, p3, p3, p3);

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		// printArrayDouble(ha0, ng, ng);
		// printArrayDouble(theta, p3, p3);
		// printArrayDouble1(hessz, p3, p3, p3);

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk;
		}
		// if(step<2){
		// 	printf("nlogelr = %f,   ologelr = %f\n", nlogelr, ologelr);
		// 	printArrayDouble1(w, n, 5, p);
		// }

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}
	// printf("step = %d\n", step);

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(xy);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

void _EstLinearCP_Multi(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *xy, *zy, *ay, *resids, *resids1, *xb;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			

			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while(step < maxstep){
		step++;

		// calculate alpha and beta
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			ologelr += yk*yk;
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k];
				}
			}			
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}


		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j] = z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk;
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk;
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(xy);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

void _EstLogisticCP_Multi(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk, mu;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *zy, *ay, *resids, *resids1, *xb;
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// printArrayDouble1(w, n, 5, p);
	// printf("p3 = %d\n",p3);

	while (step < maxstep){
		step++;

		EstLogisticCP(beta0, w, y, hess, n, p, maxstep, eps);

		// printArrayDouble(beta0, p, p);
		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			mu = 0.0;
			for(j=0;j<p1;j++){
				mu += tx[j*n+i]*beta0[j];
			}
			xb[i] 	= mu;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				mu += tmp*sh[k*n+i];
			}
			wk 	= 1.0/(1.0+exp(-mu));
			ologelr += y[i]*mu - ( mu - log(wk) );
			yk 	= y[i] - wk;
			resids[i] = yk;
			
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}
			wk = wk*(1.0-wk);
			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*wk;
				}
			}			
		}
		ologelr = -2*ologelr;

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}
		// printArrayDouble(ha, ng, ng);

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				mu 	= xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					mu += xb[n+k*n+i]*tmp;
				}
				wk 	= 1.0/(1.0+exp(-mu));
				resids1[i] = y[i] - wk;
				nlogelr += y[i]*mu - ( mu - log(wk) );
			}
			nlogelr = -2*nlogelr;
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			mu = 0.0;
			for(k=0; k < ng; k++){
				mu += xb[n+k*n+i]*dsh[k*n+i];
			}
			// wk 	= 1.0/(1.0+exp(-(mu + xb[i])));
			wk = y[i] - yk;
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*mu;
				zy[j] 	+= wz[j]*yk;				
			}
			wk *= 1-wk;
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*wk;
				}
			}
			// if(i<5){
			// 	// printArrayDouble(wz, p3, p3);
			// 	// printf("%f\n", wk*(1-wk));
			// 	printf("%f\n", yk);
			// }
		}
		// printArrayDouble1(hessz, p3, p3, p3);

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		// printArrayDouble(theta, p3, p3);
		// printArrayDouble1(hessz, p3, p3, p3);

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		// printArrayDouble(theta, p3, p3);

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			mu	=  xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				mu += xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			wk 	= 1.0/(1.0+exp(-mu));
			resids1[i] = y[i] - wk;
			nlogelr += y[i]*mu - ( mu - log(wk) );
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

void _EstPoissCP_Multi(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk, mu;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *zy, *ay, *resids, *resids1, *xb;
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// printArrayDouble1(w, n, 5, p);
	// printf("p3 = %d\n",p3);

	while (step < maxstep){
		step++;

		EstPoissCP(beta0, w, y, hess, n, p, maxstep, eps);

		// printArrayDouble(beta0, p, p);
		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			mu = 0.0;
			for(j=0;j<p1;j++){
				mu += tx[j*n+i]*beta0[j];
			}
			xb[i] 	= mu;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				mu += tmp*sh[k*n+i];
			}
			wk = exp(mu);
			ologelr += y[i]*mu - wk;
			yk 	= y[i] - wk;
			resids[i] = yk;
			
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}
			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*wk;
				}
			}			
		}
		ologelr = -2*ologelr;

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}
		// printArrayDouble(ha, ng, ng);

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				mu 	= xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					mu += xb[n+k*n+i]*tmp;
				}
				wk = exp(mu);
				resids1[i] = y[i] - wk;
				nlogelr += y[i]*mu - wk;
			}
			nlogelr = -2*nlogelr;

			// printf("nlogelr = %f, ologelr = %f\n", nlogelr, ologelr);
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			mu = 0.0;
			for(k=0; k < ng; k++){
				mu += xb[n+k*n+i]*dsh[k*n+i];
			}
			wk = y[i] - yk;
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*mu;
				zy[j] 	+= wz[j]*yk;				
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*wk;
				}
			}

			// if(i<5){
			// 	printArrayDouble(wz, p3, p3);
			// 	// printf("%f\n", wk*(1-wk));
			// 	// printf("%f\n", yk);
			// }

		}
		// printArrayDouble1(wz, n, 5, p);
		// printArrayDouble1(hessz, p3, p3, p3);

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		// printArrayDouble(theta, p3, p3);

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			mu	=  xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				mu += xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			wk = exp(mu);
			resids1[i] = y[i] - wk;
			nlogelr += y[i]*mu - wk;
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
			// break;
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

SEXP _EST_LINEAR_MULTI(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLinearCP_Multi(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_LOGISTICR_MULTI(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLogisticCP_Multi(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_POISSON_MULTI(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstPoissCP_Multi(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}


//-------------- ESTIMATION of GLM regression with multiple change-planes by bootstrap --------------------------------
void _EstLinearCP_Multi_Boot(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *xy, *zy, *ay, *resids, *resids1, *xb;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			

			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp*g[i];
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while(step < maxstep){
		step++;

		// calculate alpha and beta
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			ologelr += yk*yk*g[i];
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk*g[i];
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*g[i];
				}
			}			
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk*g[i];
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}


		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j] = z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk*g[i];
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp*g[i];
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk*g[i];
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(xy);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

void _EstLogisticCP_Multi_Boot(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk, mu;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *zy, *ay, *resids, *resids1, *xb;
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while (step < maxstep){
		step++;

		EstLogisticCP_Weight(beta0, w, y, g, hess, n, p, maxstep, eps);

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			mu = 0.0;
			for(j=0;j<p1;j++){
				mu += tx[j*n+i]*beta0[j];
			}
			xb[i] 	= mu;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				mu += tmp*sh[k*n+i];
			}
			wk 	= 1.0/(1.0+exp(-mu));
			ologelr += (y[i]*mu - ( mu - log(wk) ))*g[i];
			yk 	= y[i] - wk;
			resids[i] = yk;
			
			yk *= g[i];
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}
			wk = wk*(1.0-wk)*g[i];
			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*wk;
				}
			}			
		}
		ologelr = -2*ologelr;

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				mu 	= xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					mu += xb[n+k*n+i]*tmp;
				}
				wk 	= 1.0/(1.0+exp(-mu));
				resids1[i] = y[i] - wk;
				nlogelr += (y[i]*mu - ( mu - log(wk) ))*g[i];
			}
			nlogelr = -2*nlogelr;
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			mu = 0.0;
			for(k=0; k < ng; k++){
				mu += xb[n+k*n+i]*dsh[k*n+i];
			}
			wk = y[i] - yk;
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*mu;
				zy[j] 	+= wz[j]*yk*g[i];				
			}
			wk *= (1-wk)*g[i];
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*wk;
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}


		nlogelr = 0.0;
		for(i=0;i<n;i++){
			mu	=  xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				mu += xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			wk 	= 1.0/(1.0+exp(-mu));
			resids1[i] = y[i] - wk;
			nlogelr += (y[i]*mu - ( mu - log(wk) ))*g[i];
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

void _EstPoissCP_Multi_Boot(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk, mu;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *zy, *ay, *resids, *resids1, *xb;
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while (step < maxstep){
		step++;

		EstPoissCP_Weight(beta0, w, y, g, hess, n, p, maxstep, eps);

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			mu = 0.0;
			for(j=0;j<p1;j++){
				mu += tx[j*n+i]*beta0[j];
			}
			xb[i] 	= mu;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				mu += tmp*sh[k*n+i];
			}
			wk = exp(mu);
			ologelr += (y[i]*mu - wk)*g[i];
			yk 	= y[i] - wk;
			resids[i] = yk;
			
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk*g[i];
			}
			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*wk*g[i];
				}
			}			
		}
		ologelr = -2*ologelr;

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				mu 	= xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					mu += xb[n+k*n+i]*tmp;
				}
				wk = exp(mu);
				resids1[i] = y[i] - wk;
				nlogelr += (y[i]*mu - wk)*g[i];
			}
			nlogelr = -2*nlogelr;

			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			mu = 0.0;
			for(k=0; k < ng; k++){
				mu += xb[n+k*n+i]*dsh[k*n+i];
			}
			wk = (y[i] - yk)*g[i];
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*mu;
				zy[j] 	+= wz[j]*yk*g[i];				
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*wk;
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			mu	=  xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				mu += xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			wk = exp(mu);
			resids1[i] = y[i] - wk;
			nlogelr += (y[i]*mu - wk)*g[i];
		}
		nlogelr = -2*nlogelr;
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

SEXP _EST_LINEAR_MULTI_BOOT(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLinearCP_Multi_Boot(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_LOGISTICR_MULTI_BOOT(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLinearCP_Multi_Boot(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_POISSON_MULTI_BOOT(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth 	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];

	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstPoissCP_Multi_Boot(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

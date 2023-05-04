#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "_WAST_HEAD.h"

//-------------- Exponential mixture --------------------------------
double SST_Exp_mixture_bstr(double *y, double* thetahat, double alpha,
		int n, int K, int M, double *G, double* theta){
	int i,j,k, maxk=0, count=0;// approx=8;
	double Tn=0.0, Tn0=0.0, Tn_star, pvals, thetalpha, thetalpha1;
	double *score1, *psi, *psi_h, *Tns, psin, Vh;

	psi  	= (double*)malloc(sizeof(double)*n);
	psi_h  	= (double*)malloc(sizeof(double)*n);
	score1 	= (double*)malloc(sizeof(double)*n);
	Tns		= (double*)malloc(sizeof(double)*M);

	// if(n>500) approx = 4;
	for(i=0; i<n; i++){
		score1[i] 	= alpha*y[i] - 1.0;
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	for(k=0; k<K; k++){
		thetalpha 	= alpha/theta[k];
		thetalpha1 	= alpha - theta[k];
		for(i=0; i<n; i++){
			psi[i] = exp( thetalpha1*y[i] )/thetalpha - 1.0;
			// psi[i] = exp_approx( thetalpha1*y[i] , approx)/thetalpha - 1.0;
		}

		psin = 0.0;
		for(i=0; i<n; i++){
			psin += psi[i];
			psi_h[i] = psi[i] - score1[i]*(thetalpha - 1.0);

		}

		Vh 	= alpha*theta[k]*theta[k]*(2*theta[k] - alpha)*pow(thetalpha1, -4.0);
		Tn 	= psin*psin*Vh;

		if(Tn>Tn0){
			Tn0 = Tn;
			maxk = k;
		}
		for(j=0; j< M; j++){
			psin = 0.0;
			for(i=0; i<n; i++){
				psin += psi_h[i]*G[j*n+i];
			}
			Tn_star = psin*psin*Vh;

			if(Tn_star>Tns[j]){
				Tns[j] = Tn_star;
			}
		}

	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	thetahat[0] = theta[maxk];

	free(psi);
	free(score1);
	free(psi_h);
	free(Tns);

	return pvals;
}

SEXP _SST_Exp_mixture_bstr(SEXP Y, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, K, M;
	double alpha;

	n 		= INTEGER(DIMs)[0];
	K 		= INTEGER(DIMs)[1];
	M 		= INTEGER(DIMs)[2];

	alpha	= REAL(PARAMs)[0];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	REAL(rpvals)[0] = SST_Exp_mixture_bstr(REAL(Y), REAL(rthetahat), alpha, n, K, M, REAL(G), REAL(THETA));

	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- Normal mixture --------------------------------
double SST_Normal_mixture_bstr(double *y, double* thetahat, double alpha,
		int n, int K, int M, double *G, double* theta){
	int i,j,k, maxk=0, count=0;// approx=8;
	double tmp, Tn=0.0, Tn0=0.0, Tn_star, pvals, thetalpha;
	double *score1, *psi, *psi_h, *Tns, psin, Vh;

	psi  	= (double*)malloc(sizeof(double)*n);
	psi_h  	= (double*)malloc(sizeof(double)*n);
	score1 	= (double*)malloc(sizeof(double)*n);
	Tns		= (double*)malloc(sizeof(double)*M);

	// if(n>500) approx = 4;
	for(i=0; i<n; i++){
		score1[i] 	= alpha - y[i];
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	for(k=0; k<K; k++){
		thetalpha 	= alpha - theta[k];
		for(i=0; i<n; i++){
			psi[i] = exp( 0.5*(y[i]-alpha)*(y[i]-alpha) - 0.5*(y[i]-theta[k])*(y[i]-theta[k]) ) - 1.0;
			// psi[i] = exp_approx( 0.5*(y[i]-alpha)*(y[i]-alpha) - 0.5*(y[i]-theta[k])*(y[i]-theta[k]), approx ) - 1.0;;
		}

		psin = 0.0;
		for(i=0; i<n; i++){
			psin += psi[i];
			psi_h[i] = psi[i] - score1[i]*thetalpha;

		}

		tmp = thetalpha*thetalpha;
		Vh 	= 1.0 / ( exp(tmp) - 1.0 - tmp );
		Tn 	= psin*psin*Vh;

		if(Tn>Tn0){
			Tn0 = Tn;
			maxk = k;
		}
		for(j=0; j< M; j++){
			psin = 0.0;
			for(i=0; i<n; i++){
				psin += psi_h[i]*G[j*n+i];
			}
			Tn_star = psin*psin*Vh;

			if(Tn_star>Tns[j]){
				Tns[j] = Tn_star;
			}
		}

	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	thetahat[0] = theta[maxk];

	free(psi);
	free(score1);
	free(psi_h);
	free(Tns);

	return pvals;
}

SEXP _SST_Normal_mixture_bstr(SEXP Y, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, K, M;
	double alpha;

	n 		= INTEGER(DIMs)[0];
	K 		= INTEGER(DIMs)[1];
	M 		= INTEGER(DIMs)[2];

	alpha	= REAL(PARAMs)[0];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	REAL(rpvals)[0] = SST_Normal_mixture_bstr(REAL(Y), REAL(rthetahat), alpha, n, K, M, REAL(G), REAL(THETA));

	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- Poisson mixture --------------------------------
double SST_Poisson_mixture_bstr(double *y, double* thetahat, double alpha,
		int n, int K, int M, double *G, double* theta){
	int i,j,k, maxk=0, count=0;// approx=4;
	double tmp, Tn=0.0, Tn0=0.0, Tn_star, pvals, exptheta, thetalpha;
	double *score1, *psi, *psi_h, *Tns, psin, Vh;

	psi  	= (double*)malloc(sizeof(double)*n);
	psi_h  	= (double*)malloc(sizeof(double)*n);
	score1 	= (double*)malloc(sizeof(double)*n);
	Tns		= (double*)malloc(sizeof(double)*M);

	// if(n>500) approx = 4;
	for(i=0; i<n; i++){
		score1[i] 	= alpha - y[i];
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	for(k=0; k<K; k++){
		// exptheta 	= exp_approx(alpha-theta[k], approx);
		exptheta 	= exp(alpha-theta[k]);
		thetalpha 	= theta[k]/alpha;
		for(i=0; i<n; i++){
			psi[i] = pow(thetalpha, y[i])*exptheta - 1.0;;
		}

		psin = 0.0;
		for(i=0; i<n; i++){
			psin += psi[i];
			psi_h[i] = psi[i] - score1[i]*(1.0 - thetalpha);

		}

		tmp = alpha-theta[k];
		tmp = tmp*tmp/alpha;
		// Vh 	= 1.0 / ( exp_approx(tmp, approx) - 1.0 - tmp );
		Vh 	= 1.0 / ( exp(tmp) - 1.0 - tmp );
		Tn 	= psin*psin*Vh;

		if(Tn>Tn0){
			Tn0 = Tn;
			maxk = k;
		}
		for(j=0; j< M; j++){
			psin = 0.0;
			for(i=0; i<n; i++){
				psin += psi_h[i]*G[j*n+i];
			}
			Tn_star = psin*psin*Vh;

			if(Tn_star>Tns[j]){
				Tns[j] = Tn_star;
			}
		}

	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	thetahat[0] = theta[maxk];

	free(psi);
	free(score1);
	free(psi_h);
	free(Tns);

	return pvals;
}

SEXP _SST_Poisson_mixture_bstr(SEXP Y, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, K, M;
	double alpha;

	n 		= INTEGER(DIMs)[0];
	K 		= INTEGER(DIMs)[1];
	M 		= INTEGER(DIMs)[2];

	alpha	= REAL(PARAMs)[0];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	REAL(rpvals)[0] = SST_Poisson_mixture_bstr(REAL(Y), REAL(rthetahat), alpha, n, K, M, REAL(G), REAL(THETA));

	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- Quantile regression --------------------------------
double SST_Quantile_bstr(double *y, double *tx1, double *x, double *z, double *alphahat, double* thetahat,
		int n, int p11, int p2, int p3, double tau, int maxIter, double tol, double h, int type, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *resid1, *score1, *psi, *psin;
	double *C11, *psi_h, *K1, *Vh, *Tns, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);
	weight 	= (double*)malloc(sizeof(double)*n);


	EstQR2(tx1, y, alphahat, resid1, tau, n, p11, maxIter, tol);
	Kernelh(resid1, weight, n, h, type);

	for(i=0; i<n; i++){
		resid1[i] = (IDEX(resid1[i], 0) - tau);
		tmp 	= weight[i];
		for(j=0; j<p11; j++){
			score1[j*n+i] 	= tx1[j*n+i]*tmp;
		}
	}

	for (s = 0; s < p11; s++){
		for (t = 0; t < p11; t++){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += tx1[s*n+i]*score1[t*n+i];
			}
			C11[s*p11+t] 	= tmp/n;
		}
	}
	if(p11==1){
		C11[0] = 1.0/C11[0];
	}
	else{
		MatrixInvSymmetric(C11, p11);
	}

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
			for (s = 0; s < p2; s++){
				for (t = 0; t < p11; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						if(subg[i]){
							tmp += x[s*n+i]*tx1[t*n+i]*weight[i];
						}
					}
					K1[s*p11+t] = tmp/n;
				}
			}

			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi[s*n+i]*subg[i];
					tmp1 = 0.0;
					for (t = 0; t < p11; t++){
						tmp1 += K1[s*p11+t]*score1[t*n+i];
					}
					psi_h[s*n+i] = psi[s*n+i]*subg[i] - tmp1;

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
				Tn_star = 0.0;
				for (s = 0; s < p2; s++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi_h[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp;
				}


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
	free(resid1);
	free(score1);
	free(C11);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(Tns);
	free(weight);
	return pvals;
}

SEXP _STT_Quantile_bstr(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, K, M, maxIter, type;
	double tol, h, tau;

	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];
	type 	= INTEGER(DIMs)[7];

	tau		= REAL(PARAMs)[0];
	h 		= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rpvals, rAlpha, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));


	REAL(rpvals)[0] = SST_Quantile_bstr(REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha), REAL(rthetahat),
					n, p1, p2, p3, tau, maxIter, tol, h, type, K, M, REAL(G), REAL(THETA));
	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,  mkChar("alpha"));
	SET_STRING_ELT(list_names, 	2,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	SET_VECTOR_ELT(list, 		2, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

//-------------- Probit regression --------------------------------
void Est_probit12(int *y, double *x, double *alpha, double *residual, double *weight, double *hess, int n, int p, int maxIter, double tol){
	int i,j, k,step=0;
	double tmp, tmp1, phix1, phix2, bnorm, bnorm0 = 1.0, mui;
	double *alpha0, *dpsi, *psi10;

	alpha0 	= (double*)malloc(sizeof(double)*p);
	psi10 	= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0; j < p; j++){
		alpha0[j] = 0.0;
	}
	while(step<maxIter){
		step++;
		for(j=0; j < p; j++){
			psi10[j] = 0.0;
		}
		for(i=0; i<n; i++){
			mui = 0.0;
			for(j=0; j < p; j++){
				mui += x[j*n+i]*alpha0[j];
			}
			tmp 	= exp(-0.5*mui*mui)*MPI2;
			tmp1 	= 0.5*erf(SQRT2*mui) +0.5;
			phix1 	= tmp/tmp1;
			phix2 	= tmp/(1-tmp1);
			tmp 	= y[i]?phix1:(-phix2);
			tmp1 	= y[i]?(-phix1*(mui+phix1)):(phix2*(mui - phix2));
			for(j=0; j < p; j++){
				dpsi[j*n+i] = tmp1*x[j*n+i];
				psi10[j] 	+= tmp*x[j*n+i];
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

	for(i=0; i<n; i++){
		mui = 0.0;
		for(j=0; j < p; j++){
			mui += x[j*n+i]*alpha0[j];
		}
		tmp 	= exp(-0.5*mui*mui)*MPI2;
		tmp1 	= 0.5*erf(SQRT2*mui) +0.5;
		phix1 	= tmp/tmp1;
		phix2 	= tmp/(1-tmp1);
		tmp 	= y[i]?phix1:(-phix2);
		tmp1 	= y[i]?(-phix1*(mui+phix1)):(phix2*(mui - phix2));
		weight[i] 	= tmp1;
		residual[i] = tmp;
		for(j=0; j < p; j++){
			dpsi[j*n+i] = tmp1*x[j*n+i];
		}
	}

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*dpsi[k*n+i];
			}
			hess[j*p+k] = tmp/n;
		}
	}
	MatrixInvSymmetric(hess,p);

	free(alpha0);
	free(psi10);
	free(dpsi);
}

double SST_Probit_bstr(int *y, double *tx1, double *x, double *z, double *alphahat, double* thetahat,
		int n, int p11, int p2, int p3, int maxIter, double tol, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *resid1, *score1, *psi, *psin;
	double *C11, *psi_h, *K1, *Vh, *Tns, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);
	weight 	= (double*)malloc(sizeof(double)*n);

	Est_probit12(y, tx1, alphahat, resid1, weight, C11, n, p11, maxIter, tol);

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
		resid1[i] = resid1[i]*weight[i];
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
			for (s = 0; s < p2; s++){
				for (t = 0; t < p11; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						if(subg[i]){
							tmp += x[s*n+i]*tx1[t*n+i]*resid1[i];
						}
					}
					K1[s*p11+t] = tmp/n;
				}
			}

			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi[s*n+i]*subg[i];
					tmp1 = 0.0;
					for (t = 0; t < p11; t++){
						tmp1 += K1[s*p11+t]*score1[t*n+i];
					}
					psi_h[s*n+i] = psi[s*n+i]*subg[i] - tmp1;

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
				Tn_star = 0.0;
				for (s = 0; s < p2; s++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi_h[s*n+i]*G[j*n+i];
					}
					psin[s] = tmp;
				}


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
	free(resid1);
	free(score1);
	free(C11);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(Tns);
	free(weight);
	return pvals;
}

SEXP _STT_Probit_bstr(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, K, M, maxIter;
	double tol;

	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];
	tol 	= REAL(PARAMs)[0];

	SEXP rpvals, rAlpha, rthetahat, list, list_names;
	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));


	REAL(rpvals)[0] = SST_Probit_bstr(INTEGER(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha), REAL(rthetahat),
					n, p1, p2, p3, maxIter, tol, K, M, REAL(G), REAL(THETA));
	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,  mkChar("alpha"));
	SET_STRING_ELT(list_names, 	2,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	SET_VECTOR_ELT(list, 		2, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

//-------------- Double Robust regression --------------------------------
void EstLogisticR2(double *residual, const double *x, const double *y, double *beta, double *hess, int n, int p, int maxstep, double eps, double *weight){
	int i,j,k, step=0;
	double *beta0, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
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
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
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
		tmp = wk*(1.0-wk);
		weight[i] 	= tmp;
		residual[i] = yk;
		for(j=0;j<p;j++){
			dpsi[j*n+i] = x[j*n+i]*tmp;
		}
	}

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*dpsi[k*n+i];
			}
			hess[j*p+k] = tmp/n;
		}
	}
	MatrixInvSymmetric(hess,p);

	free(beta0);
	free(qy);
	free(dpsi);
}

void EstLinearR2(double *residual, const double *x, const double *y, double *beta, double *hess, int n, int p){
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

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			hess[j*p+k] *= n;
		}
	}

    free(xy);
}

double SST_DoubleRobust_bstr(double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *alphahat1, double *alphahat2, double* thetahat,
		int n, int p11, int p12, int p2, int p3, int maxIter, double tol, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, tmp2, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *resid1, *resid2, *score1, *score2, *psi, *psin;
	double *C11, *C22, *psi_h, *K1, *K2, *Vh, *Tns, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	score2 	= (double*)malloc(sizeof(double)*n*p12);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	K2		= (double*)malloc(sizeof(double)*p2*p12);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);
	weight 	= (double*)malloc(sizeof(double)*n);


	// printf("p11 = %d  \n ", p11);
	// printf("p2 = %d  \n ", p2);
	// printf("K = %d\n", K);
	// printf("M = %d\n", M);
	// printArrayDouble1(tx1, n, 5, p11);
	EstLinearR2(resid1, tx1, y, alphahat1, C11, n, p11);
	// printArrayDouble(alphahat1, p11, p11);
	EstLogisticR2(resid2, tx2, a, alphahat2, C22, n, p12, maxIter, tol, weight);
	// printArrayDouble(alphahat2, p12, p12);

	for(i=0; i<n; i++){
		for(j=0; j<p11; j++){
			tmp = 0.0;
			for (s = 0; s < p11; s++){
				tmp += C11[j*p11+s]*tx1[s*n+i];
			}
			score1[j*n+i] = tmp*resid1[i];
		}
		for(j=0; j<p12; j++){
			tmp = 0.0;
			for (s = 0; s < p12; s++){
				tmp += C22[j*p12+s]*tx2[s*n+i];
			}
			score2[j*n+i] = tmp*resid2[i];
		}

		for(j=0; j<p2; j++){
			psi[j*n+i] 	= x[j*n+i]*resid1[i]*resid2[i];
		}
		resid1[i] = resid1[i]*weight[i];
	}
	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	// tmp = 0.0;
	// for (j = 0; j < p3; j++){
	// 	tmp += theta[j]*theta[j];
	// }
	// printf("theta = %f\n", tmp);
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
		// printArrayDoubleInt(subg, 10, 10);
		// printf("sum of subg = %d \n", sumsb);
		if (sumsb>0){
			for (s = 0; s < p2; s++){
				for (t = 0; t < p11; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						if(subg[i]){
							tmp += x[s*n+i]*tx1[t*n+i]*resid2[i];
						}
					}
					K1[s*p11+t] = tmp/n;
				}

				for (t = 0; t < p12; t++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						if(subg[i]){
							tmp += x[s*n+i]*tx2[t*n+i]*resid1[i];
						}
					}
					K2[s*p12+t] = tmp/n;
				}
			}
			// printArrayDouble(K1, p11, p11);
			// printArrayDouble(K2, p12, p12);

			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi[s*n+i]*subg[i];
					tmp1 = 0.0;
					for (t = 0; t < p11; t++){
						tmp1 += K1[s*p11+t]*score1[t*n+i];
					}

					tmp2 = 0.0;
					for (t = 0; t < p12; t++){
						tmp2 += K2[s*p12+t]*score2[t*n+i];
					}

					psi_h[s*n+i] = psi[s*n+i]*subg[i] - tmp1 - tmp2;

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
			// Tn /= n;
			// printf("psin = %f    ", psin[0]);
			// printf("Vh = %f    ", Vh[0]);
			// printf("Tn = %f\n", Tn);
			if(Tn>Tn0){
				Tn0 = Tn;
				maxk = k;
			}
			// printf("maxk = %d\n", k);
			for(j=0; j< M; j++){
				Tn_star = 0.0;
				for (s = 0; s < p2; s++){
					tmp = 0.0;
					for(i=0; i<n; i++){
						tmp += psi_h[s*n+i]*G[j*n+i];
						// tmp += psi_h[s*n+i]*G[(j*p2+s)*n+i];
					}
					psin[s] = tmp;
				}


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
	// printf("count = %d\n", count);
	// printf("Tn = %f\n", Tn0);
	// printArrayDouble(Tns, 100, 10);

	free(subg);
	free(psi);
	free(resid1);
	free(resid2);
	free(score1);
	free(score2);
	free(C11);
	free(C22);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(K2);
	free(Tns);
	free(weight);
	return pvals;
}

SEXP _SST_DoubleRobust_bstr(SEXP Y, SEXP A, SEXP tX1, SEXP tX2, SEXP X, SEXP Z, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p11, p12, p2, p3, maxIter, K, M;
	double tol;
	n 		= INTEGER(DIMs)[0];
	p11		= INTEGER(DIMs)[1];
	p12		= INTEGER(DIMs)[2];
	p2 		= INTEGER(DIMs)[3];
	p3 		= INTEGER(DIMs)[4];
	K 		= INTEGER(DIMs)[5];
	M 		= INTEGER(DIMs)[6];
	maxIter = INTEGER(DIMs)[7];
	tol 	= REAL(PARAMs)[0];

	SEXP rpvals, rthetahat, rAlpha1, rAlpha2, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(rAlpha1		= allocVector(REALSXP, 	p11));
	PROTECT(rAlpha2		= allocVector(REALSXP, 	p12));
	PROTECT(list_names 	= allocVector(STRSXP, 	4));
	PROTECT(list 		= allocVector(VECSXP, 	4));


	REAL(rpvals)[0] = SST_DoubleRobust_bstr(REAL(Y), REAL(A), REAL(tX1), REAL(tX2),
						REAL(X), REAL(Z), REAL(rAlpha1), REAL(rAlpha2), REAL(rthetahat),
						n, p11, p12, p2, p3, maxIter, tol, K, M, REAL(G), REAL(THETA));



	SET_STRING_ELT(list_names, 	0,	mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,  mkChar("alpha1"));
	SET_STRING_ELT(list_names, 	2,  mkChar("alpha2"));
	SET_STRING_ELT(list_names, 	3,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rAlpha1);
	SET_VECTOR_ELT(list, 		2, 	rAlpha2);
	SET_VECTOR_ELT(list, 		3, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(6);
	return list;
}

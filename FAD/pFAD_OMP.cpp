#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include <iostream>
#include <mex.h>
#include <string.h>
#include <omp.h>

using namespace std;

/*compute the primal error and residual*/
__inline bool stopPri(double **X, double *sol, int m, int n, double absTol, double relTol, int display);

/*compute the dual error and resudual*/
__inline bool stopDual(double *sol, double *sol_o, double *U1, double *U2, double *U3, int m, int n, 
				double rho, double absTol, double relTol, int display);

/*compute the root of the quartic function*/
__inline double NewtonRoot(const double a, const double b, const double c, const double d);

/*get the function value*/
__inline double get_funval(const double *sol, const double *Y, const double lam, const int imgHeight, const int imgWidth);

/*method for isotropic TV*/
int tvl2_iso(double* sol, double* Y, const int imgHeight, const int imgWidth, const double lam, const double rho, 
			  const int maxIter, const double absTol, const double relTol, const int display);

			  
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* inputImg =          mxGetPr(prhs[0]);
    int imgHeight =                 (int) mxGetM(prhs[0]);
    int imgWidth =                  (int) mxGetN(prhs[0]);
    double lam =               mxGetScalar(prhs[1]);    
    double gamma =             mxGetScalar(prhs[2]);    
    int maxIter =          (int) mxGetScalar(prhs[3]);    
    double* tol =              mxGetPr(prhs[4]);    
    int display =              (int) mxGetScalar(prhs[5]);
    double *outputImg, *iter, *funVal;
    plhs[0] = mxCreateDoubleMatrix(imgHeight, imgWidth, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    outputImg = mxGetPr(plhs[0]);
	funVal = mxGetPr(plhs[1]);
    iter = mxGetPr(plhs[2]);
    memcpy(outputImg, inputImg, imgHeight * imgWidth * sizeof(double));
    iter[0] = (double)tvl2_iso(outputImg, inputImg, imgHeight, imgWidth, lam, gamma, maxIter,tol[0],tol[1],display);
	funVal[0] = get_funval(outputImg,inputImg,lam,imgHeight,imgWidth);
}
			  
__inline double NewtonRoot(const double a, const double b, const double c, const double d){
	// initial x = 0.1
	double x = 0.1;
	double x2 = 0.01;
	double x3 = 0.001;
	double x4 = 0.0001;
	double y = 1e-4 + a*1e-3 + b*1e-2 + c*0.1 + d;
	double dy = 0.;
	int count = 0;

	while(y<0){
		x*=2;
		x2 = 2*x;
		x3 = 4*x;
		x4 = 8*x;
		y = x4 + a*x3 + b*x2 + c*x + d;
	}
	while(fabs(y)>1e-15){
		dy = 4*x3 + 3*a*x2 + 2*b*x + c;
		x-=y/dy;
		x2 = x*x;
		x3 = x2*x;
		x4 = x2*x2;
		y = x4 + a*x3 + b*x2 + c*x + d;
		count ++;
		if(count > 100){
			printf("Newton method fails %f\n", (float) y);
		}
	}
	
	return x;
}

bool stopPri(double **X, double *sol, int m, int n, double absTol, double relTol, int display)
{
	const double sqrt3 = sqrt(3.0);
	double rk = 0, epri = 0, ex = 0, ez = 0;
	const int imgDim = m*n;
	int i;
	int max_threads = omp_get_max_threads();
	
	omp_set_num_threads(max_threads);
	#pragma omp parallel shared (X,sol) private(i) reduction(+: rk, ex, ez) 
	{
		#pragma omp for
		for(i = 0; i<imgDim; i++)
		{
			rk += (X[0][i]-sol[i])*(X[0][i]-sol[i]) + (X[1][i]-sol[i])*(X[1][i]-sol[i]) + (X[2][i]-sol[i])*(X[2][i]-sol[i]);
			ex += X[0][i]*X[0][i] + X[1][i]*X[1][i] + X[2][i]*X[2][i];
			ez += sol[i]*sol[i];
		}
	}
	
	
	
	rk = sqrt(rk);
	ez = sqrt3*sqrt(ez);
	ex = sqrt(ex);
	epri = ex>ez?ex:ez;
	i = m>n?m:n;
	epri *= relTol;
	epri += sqrt3*i*absTol;

	if(display == 1){
		printf("PriError and Resideual is %f %f\n", (float) epri, (float) rk);
	}
	return (rk <= epri);
}


// dual error and dual residual: stop 2
bool stopDual(double *sol, double *sol_o, double *U1, double *U2, double *U3, int m, int n, 
				double rho, double absTol, double relTol, int display)
{
	double sqrt3 = sqrt(3.0);
	double sk = 0, ed = 0;
	int imgDim = m*n;
	int i;
	int max_threads = omp_get_max_threads();
	
	omp_set_num_threads(max_threads);
	#pragma omp parallel shared(sol,sol_o,U1,U2,U3,imgDim) private(i) reduction(+: sk,ed)
	{	
		#pragma omp for
		for (i=0;i<imgDim;i++)
		{
			sk += (sol[i]-sol_o[i])*(sol[i]-sol_o[i]);
			ed += U1[i]*U1[i] + U2[i]*U2[i] + U3[i]*U3[i];
		}
	}
	
	sk =rho*sqrt3*sqrt(sk);
	ed = sqrt(ed);
	ed *= relTol;
	i = m>n?m:n;
	ed += sqrt3*i*absTol;
	
	if(display == 1){
		printf("DualError and Resideual is %f %f\n", (float) ed, (float) sk);
	}
	
	return (sk <= ed);
}

double get_funval(const double *sol, const double *Y, const double lam, const int imgHeight, const int imgWidth){
	
	double obj = 0.;
	double temp = 0.;
	int i, j;
	int imgDim = imgHeight * imgWidth;
		
	for (i=0; i<imgDim;i++){
		obj += (sol[i] - Y[i])*(sol[i] - Y[i]);
	}
	obj *=0.5;

	for(i=0;i<imgHeight-1;i++){
		for(j=0;j<imgWidth-1;j++){
			temp += sqrt((sol[i+j*imgHeight]- sol[i+j*imgHeight+1])*(sol[i+j*imgHeight]- sol[i+j*imgHeight+1]) + (sol[i+j*imgHeight]- sol[i+(j+1)*imgHeight])*(sol[i+j*imgHeight]- sol[i+(j+1)*imgHeight]));
		}			
	}
				
	for(j=0;j<imgWidth-1;j++){
		temp += fabs(sol[imgHeight-1+j*imgHeight] - sol[imgHeight + imgHeight-1+j*imgHeight]);
	}

	for(i=0;i<imgHeight-1;i++){
		temp += fabs(sol[i + imgHeight*(imgWidth-1)] - sol[1+i + imgHeight*(imgWidth-1)]);
	}
	
	return obj + lam*temp;
}

int tvl2_iso(double* sol, double* Y, const int imgHeight, const int imgWidth, const double lam, const double rho, 
			  const int maxIter, const double absTol, const double relTol, const int display)
{
	double *U1, *U2, *U3, *sol_o, *T, td, w0,w1,w2, temp1, temp2;
	double **X;
	int iter, i, j;
	unsigned int *BlkInd;
	const int imgDim = imgHeight*imgWidth;
	double wt0,wt1;
	// set some frequently used constants
	const double flam = lam/rho;
	const double rhoinv = 1.0/rho;
	const double flamp2 = flam*flam;
	const double flamp2_12 = 12*flamp2;
	const double flamp2_22 = 22*flamp2;
	const double flamp2_9 = 9*flamp2;
	const double invflamp2 = 1/flamp2;
	const double flamp4 = flamp2*flamp2;
	const double flamp6 = flamp4*flamp2;
	const double sqrt2 = sqrt(2.0);
	const double sqrt3 = sqrt(3.0);
	const double aa = 8.0*flamp2;
	const double inv3 = 1/3.;
	double alpha, temp,obj, ox1,ox2,ox3,oz;

	// intilize some used matrix
	double GGTI[4] = {0, 0,
					  0, 0};
	double GTGGTI[9] = {-1.0, 1.0, 0, 
						  0, -1.0, 1.0,
						  0, -1.0 , 1.0};
	
	double bb,cc,dd, aa2, bb2,aabb;
	double t1,t2;
	
	X = (double **) malloc(3*sizeof(double *));
	X[0] = (double *)malloc(imgDim*sizeof(double));
	X[1] = (double *)malloc(imgDim*sizeof(double));
	X[2] = (double *)malloc(imgDim*sizeof(double));
	BlkInd = (unsigned int *) malloc(imgDim*sizeof(unsigned int));
	sol_o  =  (double *) malloc(imgDim*sizeof(double));
	U1  =  (double *) calloc(imgDim, sizeof(double));
	U2  =  (double *) calloc(imgDim, sizeof(double));
	U3  =  (double *) calloc(imgDim, sizeof(double));

	int max_threads = omp_get_max_threads();
	
	// initilization
	td = 1.0/(1 + 3* rho);

	
	memcpy(sol,Y,imgDim*sizeof(double));
	memcpy(sol_o,Y,imgDim*sizeof(double));
		
	omp_set_num_threads(max_threads);
	#pragma omp parallel shared(BlkInd,imgHeight,imgWidth) private(i,j)
	{
		#pragma omp for //schedule(dynamic, NUM_CHUNK)
		for(j = 0; j<imgWidth; j++)
		{
			for(i = j%3; i<imgHeight; i+=3)
			{
				BlkInd[i+j*imgHeight] = 0;
			}

			for(i = (j+1)%3; i<imgHeight; i+=3)
			{
				BlkInd[i+j*imgHeight] = 1;
			}

			for(i = (j+2)%3; i<imgHeight; i+=3)
			{
				BlkInd[i+j*imgHeight] = 2;
			}
		}
	}
	
	int k = 0;
	int ind0 = 0;
	int ind1 = 0;
	int ind2 = 0;
	for(iter = 0; iter<maxIter; iter++)
	{
		omp_set_num_threads(max_threads);
		#pragma omp parallel shared(X,sol,U1,U2,U3) private(i) 
		{
			#pragma omp for //schedule(dynamic)
			for(i=0;i<imgDim;i++)
			{
				X[0][i] = sol[i] - U1[i] * rhoinv;
				X[1][i] = sol[i] - U2[i] * rhoinv;
				X[2][i] = sol[i] - U3[i] * rhoinv;
			}
		}		
			
		omp_set_num_threads(max_threads);
		#pragma omp parallel shared(X,BlkInd) private(i,j,k,ind0,ind1,ind2,w0,w1,w2,wt0,wt1,bb,cc,dd,alpha,temp,GGTI,GTGGTI) //schedule(dynamic, NUM_CHUNK)
		{
			#pragma omp for
			for(j = 0; j<imgWidth-1; j++)
			{
				for(i = 0; i<imgHeight-1; i++)
				{
				
					k = BlkInd[i+j*imgHeight];
					ind0 = i+1+j*imgHeight;
					ind1 = i+j*imgHeight;
					ind2 = i+(j+1)*imgHeight;

					// set [w0;w1;w2];
					w0 = X[k][ind0];w1 = X[k][ind1];w2=X[k][ind2];

					// compute [wt0;wt1] = (GG^T)^-1Gw
					wt0 = (w1 + w2 - 2*w0)*inv3;
					wt1 = (2*w2 - w1 - w0)*inv3;
					
					// compute lam_max
					if (flamp2 >= wt1*wt1 + wt0*wt0)
					{
						// compute u = (w0 + w1 + w2)/3;
						X[k][ind0] = X[k][ind1] = X[k][ind2] = (w0 + w1 + w2)*inv3; 
						continue;
					}

					////////////////////////////////////////////
					// compute w\tilde, since S\Sigma is fixed and constant
					wt0 = w0 - 2*w1 + w2;
					wt1 = w2 - w0;
					wt0 *=wt0; wt0 *=0.5;
					wt1 *=wt1; wt1 *=0.5;

					// compute c0,c1,c2,c3 in the paper
					//aa = 8.0*flamp2;
					bb = flamp2*(flamp2_22 - wt0 - wt1);
					cc = 2.0*flamp4 * (flamp2_12 - wt0 - 3*wt1);
					dd = flamp6 * (flamp2_9 - wt0 - 9*wt1);
					

					// compute alpha use the fourth solution in the closed form soltuions 
					alpha = NewtonRoot(aa,bb,cc,dd);
					
					/* compute inv(GG' + alpha/lamp2*I), note that the regualrziaton paramter is flam*/
					alpha *= invflamp2;
					temp = alpha*alpha + 4.*alpha + 3;
					GGTI[0] = GGTI[3] = (alpha + 2.)/temp;
					GGTI[1] = GGTI[2] = 1./temp;

					/* I - G'inv(GG' + alpha/lamp2*I)G*/
					GTGGTI[0] = 1.0 - GGTI[0];
					GTGGTI[4] = 1.0 - GGTI[0] + GGTI[1] + GGTI[2] - GGTI[3];
					GTGGTI[8] = 1.0 - GGTI[3];
					GTGGTI[1] = GTGGTI[3] = GGTI[0] - GGTI[1];
					GTGGTI[2] = GTGGTI[6] = GGTI[1];
					GTGGTI[5] = GTGGTI[7] = GGTI[3] - GGTI[1];

					/* compute u*/
					X[k][ind0]   = (GTGGTI[0]*w0 + GTGGTI[1]*w1 + GTGGTI[2]*w2);
					X[k][ind1]     = (GTGGTI[3]*w0 + GTGGTI[4]*w1 + GTGGTI[5]*w2);
					X[k][ind2] = (GTGGTI[6]*w0 + GTGGTI[7]*w1 + GTGGTI[8]*w2);

				}
			}

			#pragma omp for
			for(i = 0; i < imgHeight-1; i++)
			{
				k = BlkInd[i + (imgWidth-1)*imgHeight];
				ind0 = i+(imgWidth-1)*imgHeight;
				ind1 = i+1+(imgWidth-1)*imgHeight;
				w0 = X[k][ind0]; w1 = X[k][ind1];
				if(w0 > w1 + 2*flam){
					X[k][ind0] = w0 - flam; X[k][ind1] = w1 + flam;}
				else if(w1 > w0 + 2*flam){
					X[k][ind0] = w0 + flam; X[k][ind1] = w1 - flam;}
				else
				{
					X[k][ind0] = X[k][ind1]= (w0 + w1)*0.5;
				}
			}

			#pragma omp for
			for(j = 0; j<imgWidth-1; j++)
			{
				k = BlkInd[imgHeight-1+j*imgHeight];
				ind0 = imgHeight-1+j*imgHeight;
				ind1 = imgHeight-1+(j+1)*imgHeight;
				w0 = X[k][ind0]; w1 = X[k][ind1];
				if(w0 > w1 + 2*flam)
				{
					X[k][ind0] = w0 - flam; X[k][ind1] = w1 + flam;
				}
				else if(w1 > w0 + 2*flam)
				{
					X[k][ind0] = w0 + flam; X[k][ind1] = w1 - flam;
				}
				else
				{
					X[k][ind0] = X[k][ind1]= (w0 + w1)*0.5;
				}
			}
		}
			
	/*compute sol, U1, U2, U3*/
		td = 1./(1+3*rho);			
		omp_set_num_threads(max_threads);
		#pragma omp parallel shared(sol, Y, U1, U2, U3, X,td) private(i)
		{
			#pragma omp for
			for(i=0; i<imgDim;i++)
			{
				sol[i] = (Y[i] + U1[i] + U2[i] + U3[i] 
							+ rho*(X[0][i] + X[1][i] + X[2][i]))*td;
				U1[i] += rho * (X[0][i] - sol[i]);
				U2[i] += rho * (X[1][i] - sol[i]);
				U3[i] += rho * (X[2][i] - sol[i]);
			}
		}

		/*determine if the stop conditions are satisfied*/
		int step_size = 50;
		if(iter<=50)
		{
			step_size = 50;
		}else if(iter<=100){
			step_size = 20;
		}else{
			step_size = 10;
		}

		if(iter%step_size == 0){
			if(!stopDual(sol,sol_o,U1,U2,U3,imgHeight,imgWidth, rho,absTol,relTol,display)){
				memcpy(sol_o, sol, imgDim*sizeof(double));
				continue;
			}		
			if(stopPri(X,sol,imgHeight,imgWidth, absTol,relTol,display)){
				break;
			}
		}	
			
		memcpy(sol_o, sol, imgDim*sizeof(double));
	}

	
	free(X[0]);
	free(X[1]);
	free(X[2]);
	free(X);
	free(BlkInd);
	free(U1);
	free(U2);
	free(U3);
	free(sol_o);
	
	return iter;
}

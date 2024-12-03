/* file OptBlock.c
| copyright (C) 2002-2004 by Robert E. Wheeler
*/

#include "wheeler.h"
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "OptBlock.h"


int		MAXN=0;
int		nColumns=0;
bool    doWholeBlock=false; /* create block within block interactions */
bool    extraBlock=false;  /* True when candidate list is in an extra block */
bool	initRows=false;    /* True when initial design rows are specified */
bool    obScaled=false;    /* When true orthogonal blocks are scaled */

/********************************************************************** */

void printMatrix(char *name, double *matrix, int rows, int cols) {
	fprintf(stderr,"\n---------\n");
	fprintf(stderr,"%s\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(stderr,"%f ", matrix[i * cols + j]);
        }
        fprintf(stderr,"\n");
    }
	fprintf(stderr,"---------\n");
}
/* PermuteB **********************************************************
|	Randomly pemutes the n integers in a[] using the Fike
|	algorithm.  See Fike, "A permutation generation method"  The Computer
|	Journal, 18-1, Feb 75, 21-22.
*/

void PermuteB(
int		*a, 
int		n
)
{
   int 		i,
			j,
			temp;

	
   //GetRNGstate();
   for (i = 1; i < n; i++) {
	  //double unif_rand = ((float) rand() / (float)(RAND_MAX));
	  double unif_rand = 0.5;
	  //printf("unif_rand: %f\n", unif_rand);
      j =(int)((double)(1+i)*unif_rand);
	  //printf("j: %d\n", j);
      temp = a[j];
      a[j] = a[i];
      a[i] = temp;
   }
   //PutRNGstate();
}

/* RotateB ****************************************************************
|  Rotates in the new row using Given's rotations without
|  square roots.  weight is the square of the scale factor for the new row.
|  For X'(X,Y)=T'D(T,U), the result (returned in matrixXY) is (T,U) with unit
|   diagonals.  These being understood, D is placed on the diagonal.
|	NOTE: matrixXY is nTerms by nColumns, with the T part stored as an upper triangular matrix
|		requiring nTerms(nTerms+1)/2 elements. Imat indexes the elements.
|	tVec is needed if vec is to be preserved.
|	If weight is -1, vec is removed from T
|	matrixXY should initially be zero unless vec is to be added or removed from existing T
*/

#define TOLROT 1.0e-50

void RotateB(
	double	*vec,
	double	*tVec,
	double	*matrixXY,
	int		nTerms,
	int		nColumns,
	double	weight
)
{
	double	d,
			dp,
			c,
			s,
			x,
			r;
	int		i,
			j,
			kIndex;
	bool	skip;

	printMatrix("vec", vec, nColumns, 1);	
	//printMatrix(tVec, nColumns, 1);
	printMatrix("matrixXY before", matrixXY, nTerms, nColumns);
	//fprintf(stderr,"matrixXY: %2.2f\n", matrixXY);

	for (i=0;i<nColumns;i++) {
		tVec[i]=vec[i];		
	}

	skip=false;
	for (i=0;i<nTerms;i++) {
		if(! skip) {
			x=tVec[i];
			if (x==0.0)
				continue;
				
			d=matrixXY[kIndex=Imat(i,i)];
			fprintf(stderr,"i: %d, kIndex: %d, d: %2.2f, x: %2.2f\n", i, kIndex, d, x);
			dp=d+weight*x*x;
			if (fabs(dp)<TOLROT)
				continue; 
			matrixXY[kIndex]=dp;
			fprintf(stderr,"matrixXY1[kIndex]: %2.2f\n", matrixXY[kIndex]);
			c=d/dp;
			s=weight*x/dp;
			
			if(d==0.0) {
				skip=true;
				weight=0.0;
			} 				/* to avoid 0/0 since matrixXY is originally 0*/
			else 
				weight*=c;				
			kIndex++;
			for (j=i+1;j<nColumns;j++,kIndex++)	{
				r=matrixXY[kIndex];
				matrixXY[kIndex]=s*tVec[j]+c*r;
				fprintf(stderr,"matrixXY2[kIndex]: %2.2f\n", matrixXY[kIndex]);
				tVec[j]-=x*r;
			}
		}
		else
		  break;
	}
	printMatrix("matrixXY after", matrixXY, nTerms, nColumns);

}

/* BacksolveB **************************************************************
|  Backsolves Th=y where (T,y) is stored in matrixXY with T
|		upper triangular.  See RotateB() for a description of T.
|	Solves in place, resulting in (T^-1,b), where b is	the coefficients array.
|  The dimensions of (T,y) are nTerms by nColumns
|			y is a nTerms by nColumns-nTerms array.
*/

void BacksolveB(
	double	*matrixXY,
	int		nTerms,
    int		nColumns,
    int		doInverse
)
{
   int	col,
	    i,
		j,
		lIndex,
		kIndex;

	if (doInverse!=INVERSEONLY) {
			/* BacksolveB to get the coefficients */
		for (col=nColumns-1;col>=nTerms;col--) {
      	for(i=nTerms-1;i--;) {
				lIndex=Imat(i,col);
				kIndex=Imat(i,nTerms-1);
         	for (j=nTerms-1;j>i;j--)
            	matrixXY[lIndex]-=matrixXY[Imat(j,col)]*matrixXY[kIndex--];
 			}
		}
	}

	if (doInverse==NOINVERSE)
		return;

		/* BacksolveB for the inverse of T, which is upper triangular */
 	for (col=nTerms-1;col;col--) {
		for (j=col;j--;) {
		   lIndex=Imat(j,col);
		   kIndex=Imat(j,j+1);
			matrixXY[lIndex]*=(-1.0);
			for (i=1;i<col-j;i++)
				matrixXY[lIndex]-=matrixXY[Imat(i+j,col)]*
					matrixXY[kIndex++];
		}
	}

	for (i=0;i<nTerms;i++)
		matrixXY[Imat(i,i)]=1.0/matrixXY[Imat(i,i)];

}




/* getRangeB ***********************************************************************
|  increments max and min with values from a new vector
*/

void getRangeB(
	double *pMx,
	double *pMn,
	double *vec,
	int k
)
{
	int i;

	for (i=0;i<k;i++) {
		pMx[i]=maxm(pMx[i],vec[i]);
		pMn[i]=minm(pMn[i],vec[i]);
	}

}


/* Difference ***********************************************************************
|   Subtracts second vector from first
*/

void Difference(
	double	*vec,
	double	*first,
	double	*second,
	int		k
)
{
	int i;

	printMatrix("first", first, k, 1);
	printMatrix("second", second, k, 1);
	for (i=0;i<k;i++) {
		vec[i]=first[i]-second[i];
	}

}




/* reduceXtoT ********************************************************************
|	Reduces design rows in X to T via Given's rotations. Returns the determinant
*/

double reduceXtoT(
	double	*X,
	double	*T,
	int		*B,
	double  *blockMeans,
	int		k,
	int		nB,
	int     *blocksizes,
	double  *blockFactors,
	double  *vec,
	double  *Sc,
	bool    *singular
)
{
	double	*xri;
	int		*pB;
	double  *xmi;
	double logdet=0;
	double  *pMx=Sc;
	double  *pMn=Sc+k;
	double  maxT=-1e16;
	double  minT=1e16;
	double  *tvec=vec+k;
	double  t;
	double  r;
	int		K=(k*(k+1))/2;
	int		i;
	int		j;
	int		l;
	int		rowNo;
	double  *deli=0;
	double  *pT;

	*singular=false;

	for (i=0;i<k;i++) {
		pMx[i]=maxT;
		pMn[i]=minT;
	}

	memset((void *)T,0,K*sizeof(double));




	for (i=0;i<nB;i++) {
		pB=B+i*MAXN;				/* pointer to first row no in block i */
		xmi=blockMeans+i*k;			/* pointer to first mean in block i */
		if (doWholeBlock)
			deli=blockFactors+i*k;
		for (j=0;j<blocksizes[i];j++) {
			rowNo=*(pB++);
			xri=X+rowNo*k;           /* pointer to first element of row rowNo */
			Difference(vec,xri,xmi,k);
			if (doWholeBlock) {
				for (l=0;l<k;l++) {
					vec[l]*=deli[l];
				}
			}
			getRangeB(pMx,pMn,vec,k);
			RotateB(vec,tvec,T,k,k,1.0);  
			fprintf(stderr,"i: %d, j: %d\n", i, j);
		}
	}

	pT=T;
	for (i=0;i<k;i++) {
		r=(pMx[i]+pMn[i])/2;
		t=*pT;
		if (t<=0  || t<r*1e-10) {
			*singular=true;
			return(0);
		}
		logdet+=log(t); 
		pT+=k-i;
	}
	return (logdet);
}





/* reduceXtoTDpc ********************************************************************
|	For each block, reduces design rows in X to T via Given's rotations. Returns the 
|	product of the determinants for all blocks. T contains nB triangular matrices, one 
|	for each block.
*/

double reduceXtoTDpc(
	double	*X,
	double	*T,
	int		*B,
	double  *blockMeans,
	int		N,
	int		k,
	int		nB,
	int     *blocksizes,
	double  *blockFactors,
	double  *vec,
	double  *Sc,
	bool    *singular
)
{
	double	*xri;
	int		*pB;
	double  *pT;
	double  *pTT;
	double  *xmi;
	double	logdet=0;
	double  logCurDet=0;
	double  *pMx=Sc;
	double  *pMn=Sc+k;
	double  maxT=-1e16;
	double  minT=1e16;
	double  *tvec=vec+k;
	double  t;
	double  r;
	int		K=(k*(k+1))/2;
	int		curBlk;
	int		j;
	int		l;
	int		rowNo;
	int		bs;
	int		kt;
	double  *deli=0;

	*singular=false;

	for (curBlk=0;curBlk<nB;curBlk++) {
		bs=blocksizes[curBlk];

		kt=(bs<=k)?bs-1:k;

		for (j=0;j<kt;j++) {
			pMx[j]=maxT;
			pMn[j]=minT;
		}
		

		pB=B+curBlk*MAXN;				/* pointer to first row no in block curBlk */
		xmi=blockMeans+curBlk*k;			/* pointer to first mean in block curBlk */
		pT=T+curBlk*K;					/* pointer to T for block curBlk */
		memset((void *)pT,0,K*sizeof(double));

		if (doWholeBlock)
			deli=blockFactors+curBlk*k;
		for (j=0;j<bs;j++) {
			rowNo=*(pB++);
			xri=X+rowNo*k;           /* pointer to first element of row rowNo */
			Difference(vec,xri,xmi,kt);
			if (doWholeBlock) {
				for (l=0;l<kt;l++)
					vec[l]*=deli[l];
			}
			getRangeB(pMx,pMn,vec,kt);
			RotateB(vec,tvec,pT,kt,kt,1.0);  
		}
		logCurDet=0;
		pTT=pT;
		for (j=0;j<kt;j++) {
			r=(pMx[j]+pMn[j])/2;
			t=*pTT;
			if (t<=0  || t<r*1e-16) {
				*singular=true;
				return(0);
			}
			logCurDet+=log(t); 
			pTT+=kt-j;
		}
		logdet+=logCurDet/kt-log(bs);
	}


	return (logdet);
}

/* reduceXtoTDp ********************************************************************
|	For each block, reduces design rows in X to T via Given's rotations. Returns the 
|	product of the determinants for all blocks. T contains nB triangular matrices, one 
|	for each block.
*/

double reduceXtoTDp(
	double	*X,
	double	*T,
	int		*B,
	int		N,
	int		k,
	int		nB,
	int     *blocksizes,
	double  *blockFactors,
	double  *vec,
	double  *Sc,
	bool    *singular
)
{
	double	*xri;
	int		*pB;
	double  *pT;
	double  *pTT;
	double	logdet=0;
	double  logCurDet=0;
	double  *pMx=Sc;
	double  *pMn=Sc+k;
	double  maxT=-1e16;
	double  minT=1e16;
	double  *tvec=vec+k;
	double  t;
	double  r;
	int		K=(k*(k+1))/2;
	int		curBlk;
	int		j;
	int		l;
	int		rowNo;
	int		bs;
	int		kt;
	double  *deli=0;
	bool	allowSingular=false;

	*singular=false;

	for (curBlk=0;curBlk<nB;curBlk++) {
		bs=blocksizes[curBlk];

		kt=(bs>=k)?k:bs;

		allowSingular=(kt!=k); /* It is possible for xri to be zero. This is especially
							   serious when kt==bs. So return determinant from the 
							   non-zero diagonals of T. It can happen also when k==bs+1,
							   but this possibility is ignored.*/

		for (j=0;j<kt;j++) {
			pMx[j]=maxT;
			pMn[j]=minT;
		}
		

		pB=B+curBlk*MAXN;				/* pointer to first row no in block curBlk */
		pT=T+curBlk*K;					/* pointer to T for block curBlk */
		memset((void *)pT,0,K*sizeof(double));

		if (doWholeBlock)
			deli=blockFactors+curBlk*k;
		for (j=0;j<bs;j++) {
			rowNo=*(pB++);
			xri=X+rowNo*k;           /* pointer to first element of row rowNo */
			for (l=0;l<kt;l++)
				vec[l]=xri[l];
			if (doWholeBlock) {
				for (l=0;l<kt;l++) {
					vec[l]*=deli[l];
				}
			}
			getRangeB(pMx,pMn,vec,kt);
			RotateB(vec,tvec,pT,kt,kt,1.0);  
		}
		logCurDet=0;
		pTT=pT;
		for (j=0;j<kt;j++) {
			r=(pMx[j]+pMn[j])/2;
			t=*pTT;
			if (t<=0  || t<r*1e-16) {
				if (!allowSingular) { /* return product of non-zero diagonals */
					*singular=true;
					return(0);
				}
			}
			else {
				logCurDet+=log(t); 
			}
			pTT+=kt-j;
		}
		logdet+=logCurDet/kt-log(bs);
	}


	return (logdet);
}



/* makeTiFromTB **********************************************************************
|	Finds Ti, (upper triangular) the inverse of T, and then places Ti' (lower triangular) into Tip
|		so that the matrix mult XTi will be easy.
|	NOTE: T has scale values on diagonal which are incorporated in Ti
*/

void makeTiFromTB(
	double	*Tip,
	double	*T,
	double	*W,	
	double  *avVar,
	int		k
)
{	int	i;
	int j;
	int g;
	int	sizeT=sizeof(double)*(k*(k+1)/2);
	double d;
	double *pTip;
	double t;
	double aVar;
	

	W=(double *)memcpy((void*)W,(void*)T,sizeT);

	BacksolveB(W,k,k,INVERSEONLY);


	/* Put transpose of W in Tip */
	g=0;
	for (j=0;j<k;j++) {
		for (i=0;i<=j;i++) {
			Tip[g++]=W[Imat(i,j)];
		}
	}

		/* Scale rows by diagonal (DT)^-1=T^-1D-1 */
	memset((void *)W,0,k*sizeof(double)); /* reuse W to get average  variance */
	pTip=Tip;
	for (i=0;i<k;i++) {
		d=sqrt(*(pTip+i));
		*(pTip+i)=1.0;
		for (j=0;j<=i;j++) {
			t=d*(*pTip);
			*(pTip++)=t;
			W[j]+=t*t;
		}
	}

	aVar=0;
	for (i=0;i<k;i++) {
		aVar+=log(W[i]);
	}
	*avVar=exp(aVar/(double)k); /* average variance */

}

/* makeTiFromTDpc **********************************************************************
|	Finds Ti, (upper triangular) the inverse of T, and then places Ti' (lower triangular) into Tip
|		so that the matrix mult XTi will be easy.
|	NOTE: T has scale values on diagonal which are incorporated in Ti
|	This function does this for each block in T
*/



void makeTiFromTDpc(
	double	*Tip,
	double	*T,
	double	*W,	
	int     *blocksizes,
	int		nB,
	int		curBlock,
	int		newBlock,
	int		k
)
{	int	i;
	int j;
	int g;
	int l;
	int K=(k*(k+1))/2;
	int	sizeT=sizeof(double)*K;
	bool twoBlocks=false;
	double d;
	double *pTip;
	double *ppTip;
	double *pT;
	double *PTip;
	int		bs;
	int		kt;

	if (curBlock!=-1)
		twoBlocks=true;  /* only process the two specified blocks */

	for (l=0;l<nB;l++) {
		if (!twoBlocks || (curBlock==l || newBlock==l)) {
			pT=T+l*K;
			PTip=Tip+l*K;
			bs=blocksizes[l];
			kt=(bs<k)?bs:k;

			W=(double *)memcpy((void*)W,(void*)pT,sizeT);


			BacksolveB(W,kt,kt,INVERSEONLY);

			/* Put transpose of W in Tip */
			g=0;
			for (j=0;j<kt;j++) {
				for (i=0;i<=j;i++) {
					PTip[g++]=W[Imat(i,j)];
				}
			}

				/* Scale rows by diagonal */
			pTip=PTip;
			for (i=0;i<kt;i++) {
				d=sqrt(*(pTip+i));
				*(pTip+i)=1.0;
				for (j=0;j<=i;j++) {
					ppTip=pTip;
					*(pTip++)=d*(*ppTip);
				}
			}
		}
	}

}

/* makeTiFromTDp **********************************************************************
|	Finds Ti, (upper triangular) the inverse of T, and then places Ti' (lower triangular) into Tip
|		so that the matrix mult XTi will be easy.
|	NOTE: T has scale values on diagonal which are incorporated in Ti
|	This function does this for each block in T
*/



void makeTiFromTDp(
	double	*Tip,
	double	*T,
	double	*W,	
	int     *blocksizes,
	int		nB,
	int		curBlock,
	int		newBlock,
	int		k
)
{	int	i;
	int j;
	int g;
	int l;
	int K=(k*(k+1))/2;
	int	sizeT=sizeof(double)*K;
	bool twoBlocks=false;
	double d;
	double *pTip;
	double *ppTip;
	double *pT;
	double *PTip;
	int		bs;
	int		kt;

	if (curBlock!=-1)
		twoBlocks=true;  /* only process the two specified blocks */

	for (l=0;l<nB;l++) {
		if (!twoBlocks || (curBlock==l || newBlock==l)) {
			pT=T+l*K;
			PTip=Tip+l*K;
			bs=blocksizes[l];
			kt=(bs<k)?bs:k;

			W=(double *)memcpy((void*)W,(void*)pT,sizeT);


			BacksolveB(W,kt,kt,INVERSEONLY);

			/* Put transpose of W in Tip */
			g=0;
			for (j=0;j<kt;j++) {
				for (i=0;i<=j;i++) {
					PTip[g++]=W[Imat(i,j)];
				}
			}

				/* Scale rows by diagonal */
			pTip=PTip;
			for (i=0;i<kt;i++) {
				d=sqrt(*(pTip+i));
				*(pTip+i)=1.0;
				for (j=0;j<=i;j++) {
					ppTip=pTip;
					*(pTip++)=d*(*ppTip);
				}
			}
		}
	}

}



/* exchangeBlock *****************************************************************************
|	Exchanges rows between blocks
*/

void exchangeBlock(
	double	*T,
	double  *X,
	double  *vec,
	double	*blockMeans,
	int		*B,
	int		*blocksizes,
	int		xcur,
	int		xnew,
	int		curBlock,
	int		newBlock,
	int		nB,
	int		k
	)
{
	double	*xmi;
	double	*xmj;
	double	*xri;
	double	*xrj;
	double  *tvec=vec+k;
	int		rowNoi;
	int		rowNoj;
	int		ni;
	int		nj;
	double  C;
	int		i;
	int		iBlock=nB*MAXN;

	xmi=blockMeans+curBlock*k;
	rowNoi=B[IB(curBlock,xcur)];
	xri=X+rowNoi*k;
	ni=blocksizes[curBlock];



	if (extraBlock && newBlock==nB) {
		rowNoj=B[iBlock+xnew];
		xrj=X+rowNoj*k;

		for (i=0;i<k;i++)
			vec[i]=xrj[i]-xmi[i];
		RotateB(vec,tvec,T,k,k,1.0); 

		for (i=0;i<k;i++)
			vec[i]=xri[i]-xmi[i];
		RotateB(vec,tvec,T,k,k,-1.0); 
		
		for (i=0;i<k;i++)
			vec[i]=xrj[i]-xri[i];
		RotateB(vec,tvec,T,k,k,-1.0/(double)ni);


		B[iBlock+xnew]=rowNoi;

	}
	else {
		rowNoj=B[IB(newBlock,xnew)];
		xrj=X+rowNoj*k;
		xmj=blockMeans+newBlock*k;
		nj=blocksizes[newBlock];
		C=(double)(ni+nj)/(double)(ni*nj);


		for (i=0;i<k;i++)
			vec[i]=xmj[i]-xmi[i];

		RotateB(vec,tvec,T,k,k,1.0); 
		
		for (i=0;i<k;i++)
			vec[i]-=xrj[i]-xri[i];

		RotateB(vec,tvec,T,k,k,-1.0); 

		for (i=0;i<k;i++)
			vec[i]=xrj[i]-xri[i];

		RotateB(vec,tvec,T,k,k,1.0-C); 


		for (i=0;i<k;i++)
			xmi[i]+=(xrj[i]-xri[i])/(double)ni;
		for (i=0;i<k;i++)
			xmj[i]+=(xri[i]-xrj[i])/(double)nj;
		B[IB(newBlock,xnew)]=rowNoi;
	}

	B[IB(curBlock,xcur)]=rowNoj;

}

/* exchangeBlockWhole *****************************************************************************
|	Exchanges rows between blocks when doWholeBlock is true
*/

void exchangeBlockWhole(
	double	*T,
	double  *X,
	double  *vec,
	double	*blockMeans,
	int		*B,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		xnew,
	int		curBlock,
	int		newBlock,
	int		nB,
	int		k
	)
{
	double	*xmi=0;
	double	*xmj=0;
	double	*xri=0;
	double	*xrj=0;
	double  *tvec=vec+k;
	int		rowNoi;
	int		rowNoj;
	int		ni=0;
	int		nj=0;
	int		i;
	double  *deli;
	double  *delj;
	int     iBlock=nB*MAXN;


	xmi=blockMeans+curBlock*k;
	rowNoi=B[IB(curBlock,xcur)];
	xri=X+rowNoi*k;
	if (extraBlock && newBlock==nB) 
		rowNoj=B[iBlock+xnew];
	else
		rowNoj=B[IB(newBlock,xnew)];
	xrj=X+rowNoj*k;
	ni=blocksizes[curBlock];
	deli=blockFactors+curBlock*k;

		/* Swap in current block */
	for (i=0;i<k;i++)
		vec[i]=deli[i]*(xrj[i]-xmi[i]);
	RotateB(vec,tvec,T,k,k,1.0); 

	for (i=0;i<k;i++)
		vec[i]=deli[i]*(xri[i]-xmi[i]);
	RotateB(vec,tvec,T,k,k,-1.0); 

	for (i=0;i<k;i++)
		vec[i]=deli[i]*(xrj[i]-xri[i]);
	RotateB(vec,tvec,T,k,k,-1/(double)ni); 



	if (!extraBlock || newBlock!=nB) { /* Swap in new block unless it is the extra block */
		xmj=blockMeans+newBlock*k;
		nj=blocksizes[newBlock];
		delj=blockFactors+newBlock*k;	

		for (i=0;i<k;i++)
			vec[i]=delj[i]*(xrj[i]-xmj[i]);
		RotateB(vec,tvec,T,k,k,-1.0); 

		for (i=0;i<k;i++)
			vec[i]=delj[i]*(xri[i]-xmj[i]);
		RotateB(vec,tvec,T,k,k,1.0); 

		for (i=0;i<k;i++)
			vec[i]=delj[i]*(xrj[i]-xri[i]);
		RotateB(vec,tvec,T,k,k,-1/(double)nj); 
	}



	for (i=0;i<k;i++) /* Fix current block means */
		xmi[i]+=(xrj[i]-xri[i])/(double)ni;

	B[IB(curBlock,xcur)]=rowNoj;

	if (extraBlock && newBlock==nB) { /* the extra block is not mean centered */
		B[iBlock+xnew]=rowNoi;
	}
	else { /* Fix new block means */
		for (i=0;i<k;i++)
			xmj[i]+=(xri[i]-xrj[i])/(double)nj;
		B[IB(newBlock,xnew)]=rowNoi;
	}

}


/* exchangeOB *****************************************************************************
|	Exchanges rows between blocks
*/

void exchangeOB(
	double  *X,
	double  *vec,
	double  *blockMeans,
	double  *gMean,
	double  *gSS,
	double  *blockSS,
	int		*B,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		xnew,
	int		curBlock,
	int		newBlock,
	int		nB,
	int		k,
	int		Nxb
	)
{
	double  *xmi=0;
	double  *xmj=0;
	double	*xri;
	double	*xrj;
	double  *xmn=0;
	double  *gMeanOriginal=gMean+2*k;  /* gMean+k is blockSS */
	double  tSS;
	int		rowNoi;
	int		rowNoj;
	int		ni;
	int		nj;
	int		i;
	int		l;
	int     iBlock=nB*MAXN;
	double  *deli=0;
	double  *delj=0;
	double  curSS=0;
	double  t;


	rowNoi=B[IB(curBlock,xcur)];
	xri=X+rowNoi*k;  /* move from curBlock to newBlock */

	if (doWholeBlock)
		deli=blockFactors+curBlock*k;

	for (i=0;i<k;i++)
		gMeanOriginal[i]=gMean[i];  /* save a copy */

	if (extraBlock && newBlock==nB) { 
		rowNoj=B[iBlock+xnew];
		xrj=X+rowNoj*k;     /* move from extra block to curBlock */
		for (i=0;i<k;i++) {
			t=(xrj[i]-xri[i])/(double)Nxb;
			if (doWholeBlock)
				t*=deli[i];
			gMean[i]=gMeanOriginal[i]+t; /* adjust grand mean	 */
		}
	}
	else {
		rowNoj=B[IB(newBlock,xnew)];
		xmj=blockMeans+newBlock*k;
		xrj=X+rowNoj*k;  /* move from newBlock to curBlock */
	}

	xmi=blockMeans+curBlock*k;

	ni=blocksizes[curBlock];
	
	for (i=0;i<k;i++) {
		t=(xrj[i]-xri[i])/ni;
		if (doWholeBlock)
			t*=deli[i];
		xmi[i]+=t;         /* adjust current block mean */
	}


	B[IB(curBlock,xcur)]=rowNoj;
		
	if (!extraBlock || newBlock!=nB) {

		curSS=0;
		for (i=0;i<k;i++) {
			t=xmi[i]-gMean[i];
			if (doWholeBlock)
				t*=deli[i];
			if (obScaled)
				curSS+=t*t/gSS[i];
			else
				curSS+=t*t;

		}
		blockSS[curBlock]=curSS;

		nj=blocksizes[newBlock];

		xmj=blockMeans+newBlock*k; 
		
		if (doWholeBlock) {
			delj=blockFactors+newBlock*k;
		}


		curSS=0;
		for (i=0;i<k;i++) {
			t=(xri[i]-xrj[i])/nj;
			if (doWholeBlock)
				t*=delj[i];
			xmj[i]+=t;			/* adjust new block mean */
			if (obScaled)
				curSS+=xmj[i]*xmj[i]/gSS[i];
			else
				curSS+=xmj[i]*xmj[i];
		}

		blockSS[newBlock]=curSS;

		B[IB(newBlock,xnew)]=rowNoi;
	}
	else {
			/* redo all block SS because gMean has changed. */
		for (i=0;i<nB;i++) {
			tSS=0;
			xmn=blockMeans+i*k;
			for (l=0;l<k;l++) {
				t=xmn[l]+gMeanOriginal[l]-gMean[l];
				if (obScaled)
					tSS+=t*t/gSS[i];
				else
					tSS+=t*t;	
			}
			blockSS[i]=tSS;
		}

		B[iBlock+xnew]=rowNoi;
	}
}



/* exchangeDpc *****************************************************************************
|	Exchanges rows between blocks
*/

void exchangeDpc(
	double	*T,
	double  *X,
	double  *vec,
	double  *blockMeans,
	int		*B,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		xnew,
	int		curBlock,
	int		newBlock,
	int		nB,
	int		k
	)
{
	double  *xmi=0;
	double  *xmj=0;
	double	*xri;
	double	*xrj;
	double  *pT;
	double  *tvec=vec+k;
	int		rowNoi;
	int		rowNoj;
	int		ni;
	int		nj;
	int     K=(k*(k+1))/2;
	int     kc;
	int		kn;
	int		i;
	int     iBlock=nB*MAXN;
	double  *deli=0;


	xmi=blockMeans+curBlock*k;

	rowNoi=B[IB(curBlock,xcur)];
	xri=X+rowNoi*k;  /* move from curBlock to newBlock */
	if (extraBlock && newBlock==nB)
		rowNoj=B[iBlock+xnew];
	else {
		rowNoj=B[IB(newBlock,xnew)];
		xmj=blockMeans+newBlock*k;
	}
	xrj=X+rowNoj*k;  /* move from newBlock to curBlock */
	ni=blocksizes[curBlock];
	

	if (doWholeBlock) 
		deli=blockFactors+curBlock*k;
	

	kc=(ni<=k)?ni-1:k;

	pT=T+curBlock*K;

	Difference(vec,xrj,xmi,kc);
	if (doWholeBlock) {
		for (i=0;i<kc;i++)
			vec[i]=deli[i];
	}

	RotateB(vec,tvec,pT,kc,kc,1.0); /* add xrj */

	Difference(vec,xri,xmi,kc);
	if (doWholeBlock) {
		for (i=0;i<kc;i++)
			vec[i]=deli[i];
	}

	RotateB(vec,tvec,pT,kc,kc,-1.0); /* remove xri */

	Difference(vec,xrj,xri,kc);
	if (doWholeBlock) {
		for (i=0;i<kc;i++)
			vec[i]=deli[i];
	}

	RotateB(vec,tvec,pT,kc,kc,-1/(double)ni); /* adjust */
	B[IB(curBlock,xcur)]=rowNoj;

	for (i=0;i<kc;i++)				/* fix mean */
		xmi[i]+=(xrj[i]-xri[i])/(double)ni;
		
	if (!extraBlock || newBlock!=nB) {
		nj=blocksizes[newBlock];
		kn=(nj<=k)?nj-1:k;

		pT=T+newBlock*K;

		Difference(vec,xri,xmj,kn);
		if (doWholeBlock) {
			for (i=0;i<kn;i++)
				vec[i]=deli[i];
		}

		RotateB(vec,tvec,pT,kn,kn,1.0); /* add xri */

		Difference(vec,xrj,xmj,kn);
		if (doWholeBlock) {
			for (i=0;i<kn;i++)
				vec[i]=deli[i];
		}

		RotateB(vec,tvec,pT,kn,kn,-1.0); /* remove xrj */

		Difference(vec,xri,xrj,kn);
		if (doWholeBlock) {
			for (i=0;i<kn;i++)
				vec[i]=deli[i];
		}

		RotateB(vec,tvec,pT,kn,kn,-1/(double)nj); /* adjust */
	
		for (i=0;i<kn;i++)				/* fix mean */
			xmj[i]+=(xri[i]-xrj[i])/(double)nj;

		B[IB(newBlock,xnew)]=rowNoi;
	}
	else {
		B[iBlock+xnew]=rowNoi;
	}
}

/* exchangeDp *****************************************************************************
|	Exchanges rows between blocks
*/

void exchangeDp(
	double	*T,
	double  *X,
	double  *vec,
	int		*B,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		xnew,
	int		curBlock,
	int		newBlock,
	int		nB,
	int		k
	)
{
	double	*xri;
	double	*xrj;
	double  *pT;
	double  *tvec=vec+k;
	int		rowNoi;
	int		rowNoj;
	int		ni;
	int		nj;
	int     K=(k*(k+1))/2;
	int     kc;
	int		kn;
	int		i;
	int     iBlock=nB*MAXN;
	double  *deli=0;
	double  *delj=0;



	rowNoi=B[IB(curBlock,xcur)];
	xri=X+rowNoi*k;  /* move from curBlock to newBlock */
	if (extraBlock && newBlock==nB)
		rowNoj=B[iBlock+xnew];
	else 
		rowNoj=B[IB(newBlock,xnew)];
	xrj=X+rowNoj*k;  /* move from newBlock to curBlock */
	ni=blocksizes[curBlock];

	if (doWholeBlock) 
		deli=blockFactors+curBlock*k;
	

	kc=(ni>=k)?k:ni;

	pT=T+curBlock*K;

	for (i=0;i<kc;i++)
		vec[i]=xrj[i];
	if (doWholeBlock) {
		for (i=0;i<kc;i++)
			vec[i]=deli[i];
	}

	RotateB(vec,tvec,pT,kc,kc,1.0); /* add xrj */
	
	for (i=0;i<kc;i++)
		vec[i]=xri[i];
	if (doWholeBlock) {
		for (i=0;i<kc;i++)
			vec[i]=deli[i];
	}

	RotateB(vec,tvec,pT,kc,kc,-1.0); /* remove xri */
		
	if (!extraBlock || newBlock!=nB) {
		nj=blocksizes[newBlock];
		kn=(nj>=k)?k:nj;

		if (doWholeBlock) 
			delj=blockFactors+newBlock*k;

		pT=T+newBlock*K;

		for (i=0;i<kn;i++)
			vec[i]=xri[i];
		if (doWholeBlock) {
			for (i=0;i<kn;i++)
				vec[i]=delj[i];
		}
		RotateB(vec,tvec,pT,kn,kn,1.0); /* add xri */
		
		for (i=0;i<kn;i++)
			vec[i]=xrj[i];
		if (doWholeBlock) {
			for (i=0;i<kn;i++)
				vec[i]=delj[i];
		}

		RotateB(vec,tvec,pT,kn,kn,-1.0);  /* remove xrj */

	}
	B[IB(curBlock,xcur)]=rowNoj;
	if (extraBlock && newBlock==nB)
		B[iBlock+xnew]=rowNoi;
	else 
		B[IB(newBlock,xnew)]=rowNoi;
}

/* findDeltaBlock ************************************************************************
| Calculates the improvement (delta) in the design criterion when exchanging a point 
| from the current block with a point from another block.
| 
| Parameters:
| tX - Transformed design points
| tBlockMeans - Transformed block means
| B - Block assignments
| nB - Number of blocks
| nEx - Number of extra points
| blocksizes - Array containing size of each block
| xcur - Current point index to potentially exchange
| xnew - Will store index of best point to exchange with (output)
| curBlock - Current block index
| newBlock - Will store index of best block to exchange with (output)
| k - Number of dimensions/variables
*/ 

#define deltaTol 1e-12  // Minimum improvement threshold

double findDeltaBlock(
	double *tX,
	double *tBlockMeans, 
	int    *B,
	int    nB,
	int    nEx,
	int    *blocksizes,
	int    xcur,
	int    *xnew,
	int    curBlock,
	int    *newBlock,
	int    k
	)
{
	double delta=0;  // Tracks best improvement found
	double d;        // Current improvement being evaluated
	double Gi[3];    // Geometric coefficients
	double Mi[3];    // Moment terms
	double M1i[3];   // Combined geometric and moment terms
	double g, h;     // Temporary calculation variables
	double dif, dif1, dif2;  // Difference terms
	int    ni, nj;   // Block sizes
	double *fi;      // Pointer to current point to exchange
	double *fj;      // Pointer to candidate point to exchange with
	double *fmi;     // Mean of current block
	double *fmj;     // Mean of candidate block
	int i, j, l;
	int curRowNo;
	int rowNo;
	int iBlock=nB*MAXN;

	// Initialize geometric coefficients
	Gi[1]=1;
	Gi[2]=0;
	curRowNo=B[IB(curBlock,xcur)];
	ni=blocksizes[curBlock];

	fi=tX+curRowNo*k;
	fmi=tBlockMeans+curBlock*k;

	// Loop through all blocks except current
	for (i=0;i<nB;i++) {
		if (i!=curBlock) {
			nj=blocksizes[i];
			// Calculate geometric coefficient based on block sizes
			Gi[0]=(double)(ni+nj)/(double)(ni*nj);

			fmj=tBlockMeans+i*k;
			
			// Calculate squared distance between block means
			g=0;
			for (l=0;l<k;l++) {
				dif=fmj[l]-fmi[l];
				g+=dif*dif;
			}
			Mi[0]=g;

			// Try exchanging with each point in candidate block
			for (j=0;j<nj;j++) {
				rowNo=B[IB(i,j)];
				fj=tX+rowNo*k;
				g=0;
				h=0;
				
				// Calculate cross terms between means and points
				for (l=0;l<k;l++) {
					dif1=fmj[l]-fmi[l];
					dif2=fj[l]-fi[l];
					g+=dif1*dif2;
					h+=dif2*dif2;
				}
				Mi[1]=g;
				Mi[2]=h;
			
				// Combine geometric and moment terms
				for (l=0;l<3;l++)
					M1i[l]=Gi[l]+Mi[l];

				// Calculate improvement in criterion
				d=-(1+M1i[0]*M1i[2]-M1i[1]*M1i[1]); 
				
				// Update best exchange if improvement is large enough
				if ((d-delta)>deltaTol) {
					delta=d;
					*newBlock=i;
					*xnew=j;
				}
			}
		}
	}

	// Check exchanges with extra block if it exists
	if (extraBlock) {
		// Special geometric coefficients for extra block
		Gi[0]=(double)(ni+1)/(double)ni;
		Gi[1]=1/(double)ni;
		Gi[2]=-(double)(ni-1)/(double)ni;
		
		// Calculate squared distance to current point
		g=0;
		for (l=0;l<k;l++) {
			dif=fi[l]-fmi[l];
			g+=dif*dif;
		}
		Mi[2]=g;

		// Try exchanging with each extra point
		for (j=0;j<nEx;j++) {
			rowNo=B[iBlock+j];
			fj=tX+rowNo*k;
			g=h=0;
			for (l=0;l<k;l++) {
				dif1=fi[l]-fmi[l];
				dif2=fj[l]-fmi[l];
				g+=dif1*dif2;
				h+=dif2*dif2;
			}
			Mi[0]=h;
			Mi[1]=g;
			for (l=0;l<3;l++)
				M1i[l]=Gi[l]+Mi[l];

			// Calculate improvement for extra block exchange
			d=-(1+(M1i[0]*M1i[2]-M1i[1]*M1i[1]));
			if ((d-delta)>deltaTol) {
				delta=d;
				*newBlock=nB;
				*xnew=j;
			}
		}
	}
	return(delta);
}

/* transW ***************************************************************************
|  Transforms a column of W with Tip. 
|	Returns the sum of squares of the transformed W, and the transformed W itself
|		is returned in W.
*/

double transW(
	double *Tip,
	double *tVec,
	double *W,
	int k
) 
{
	double d=0;
	double *pTip;
	double g;
	int i;
	int j;

	pTip=Tip;
	for (i=0;i<k;i++) {
		g=0;
		for (j=0;j<=i;j++) {
			g+=W[j]*(*pTip++);
		}
		tVec[i]=g;
		d+=g*g;
	}
	memcpy((void*)W,(void*)tVec,k*sizeof(double));

	return d;
}



/* findDeltaBlockWhole ************************************************************************
| Calculates the delta for exchanging xold in the current block with xnew in a new block
|	Assumes block, within block interactions.
*/ 


double findDeltaBlockWhole(
	double *X,
	double *Tip,
	double *W,
	double *blockMeans,
	int		*B,
	int		nB,
	int		nEx,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		*xnew,
	int		curBlock,
	int		*newBlock,
	int		k
	)
{
	double delta=0;
	double d;
	double Gi[3]; /* Gi and Gj inverses */
	double Gj[3];
	double dii[3];
	double djj[3];
	double dij[4]; /* not symmetric */
	double Deltai[3];
	double Deltaj[3];
	double DeltaiInv[3];
	double DELTA[3];
	double det;
	int    ni;
	int    nj;
	double *xri; /* current point */
	double *xrj; /* new point */
	double *xmi; /* current block mean */
	double *xmj; /* new block mean */
	double *Wi1=W;		/* First column of Wi  */
	double *Wi2=Wi1+k;	/* Second column of Wi */
	double *Wj1=Wi2+k;	/* First column of Wj */
	double *Wj2=Wj1+k;	/* Second column of Wj */
	double *tVec=Wj2+k;  /* temp vector */
	double *deli;		/* Diagonal for i */
	double *delj;       /* Diagonal for j */
	int i;
	int j;
	int l;
	int b;
	int curRowNo;
	int rowNo;
	int iBlock=nB*MAXN;

	curRowNo=B[IB(curBlock,xcur)];
	ni=blocksizes[curBlock];
	Gi[0]=(double)(ni+1)/(double)ni;
	Gi[1]=1/(double)ni;
	Gi[2]=-(double)(ni-1)/(double)ni;

	xri=X+curRowNo*k;
	xmi=blockMeans+curBlock*k;
	deli=blockFactors+curBlock*k;

	for (i=0;i<k;i++) {
		Wi2[i]=deli[i]*(xri[i]-xmi[i]);
	}


		/* dii[3]=Wi2^T(M^-1)Wi2 */

	dii[2]=transW(Tip,tVec,Wi2,k);

	for (b=0;b<nB;b++) {
		if (b!=curBlock) {
			nj=blocksizes[b];
			Gj[0]=-(double)(nj-1)/(double)nj;
			Gj[1]=1/(double)nj;
			Gj[2]=(double)(nj+1)/(double)nj;

			xmj=blockMeans+b*k;
			delj=blockFactors+b*k;

			for (j=0;j<nj;j++) {
				rowNo=B[IB(b,j)];
				xrj=X+rowNo*k;

				 /* Make W */
				for (l=0;l<k;l++) {
					Wi1[l]=deli[l]*(xrj[l]-xmi[l]);
					Wj2[l]=delj[l]*(xri[l]-xmj[l]);
					Wj1[l]=delj[l]*(xrj[l]-xmj[l]);
				}
				/* transform W=(T')^-1 W */
				dii[0]=transW(Tip,tVec,Wi1,k);
				djj[2]=transW(Tip,tVec,Wj2,k);
				djj[0]=transW(Tip,tVec,Wj1,k);
				
				dii[1]=djj[1]=0;
				dij[0]=dij[1]=dij[2]=dij[3]=0;

				/* Finish up dii, djj, dij */
				for (l=0;l<k;l++) {
					dii[1]+=Wi1[l]*Wi2[l];
					djj[1]+=Wj1[l]*Wj2[l];
					dij[0]+=Wi1[l]*Wj1[l];
					dij[1]+=Wi1[l]*Wj2[l];
					dij[2]+=Wi2[l]*Wj1[l];
					dij[3]+=Wi2[l]*Wj2[l];
				}
				/* form Delta's */
				for (l=0;l<3;l++) {
					Deltai[l]=Gi[l]+dii[l];
					Deltaj[l]=Gj[l]+djj[l];
				}
				/* Get Deltai inverse */
				det=Deltai[0]*Deltai[2]-Deltai[1]*Deltai[1];
				DeltaiInv[0]=Deltai[2]/det;
				DeltaiInv[2]=Deltai[0]/det;
				DeltaiInv[1]=-Deltai[1]/det;
			
				/* form DELTA=Deltaj-dji(DeltaiInv)dij Note: dji=dij' */

				tVec[0]=dij[0]*DeltaiInv[0]+dij[2]*DeltaiInv[1];
				tVec[1]=dij[0]*DeltaiInv[1]+dij[2]*DeltaiInv[2];

				DELTA[0]=Deltaj[0]-tVec[0]*dij[0]-tVec[1]*dij[2];
				DELTA[1]=Deltaj[1]-tVec[0]*dij[1]-tVec[1]*dij[3];

				tVec[0]=dij[1]*DeltaiInv[0]+dij[3]*DeltaiInv[1];
				tVec[1]=dij[1]*DeltaiInv[1]+dij[3]*DeltaiInv[2];

				DELTA[2]=Deltaj[2]-tVec[0]*dij[1]-tVec[1]*dij[3];

				d=(det*(DELTA[0]*DELTA[2]-DELTA[1]*DELTA[1])-1);

				if ((d-delta)>deltaTol) {
					delta=d;
					*newBlock=b;
					*xnew=j;
				}

			}

		}
	}

	if (extraBlock) { /* See if extra rows can be swapped  */
		for (j=0;j<nEx;j++) {
			rowNo=B[iBlock+j];
			xrj=X+rowNo*k;
			 /* Make W */
			for (l=0;l<k;l++) 
				Wi1[l]=deli[l]*(xrj[l]-xmi[l]);
			/* transform W=(T')^-1 W */
			dii[0]=transW(Tip,tVec,Wi1,k);
			dii[1]=0;

			/* Finish up dii */
			for (l=0;l<k;l++) 
				dii[1]+=Wi1[l]*Wi2[l];
			/* form Delta */
			for (l=0;l<3;l++) 
				Deltai[l]=Gi[l]+dii[l];
		    d=(-(Deltai[0]*Deltai[2]-Deltai[1]*Deltai[1])-1);
			if ((d-delta)>deltaTol) {
				delta=d;
				*newBlock=nB;
				*xnew=j;
			}
		}

	}
	return(delta);
}

/* findDeltaOB ************************************************************************
| Finds a swap that reduces the SS for exchanging xold in the current block with xnew in a new block
*/ 

double findDeltaOB(
	double *X,
	double *blockMeans,
	double *vec,
	double *blockSS,
	int		*B,
	int		nB,
	int		nEx,
	int		*blocksizes,
	double  *blockFactors,
	double  *gMean,
	double  *gSS,
	int		xcur,
	int		*xnew,
	int		curBlock,
	int		*newBlock,
	int		k,
	int		Nxb,
	bool    *failure
	)
{
	double delta=0;
	double *xri=0;
	double *xrj=0;
	double *xmi=0; 
	double *xmj=0; 
	double *xmn=0;
	double *xrc=vec+k;
	double *gMeanSwap=gMean+2*k; /* blockSS is at gMean+k */

	int newBlk;
	int xnw=0;
	int i;
	int j;
	int l;
	int curRowNo;
	int newRowNo;
	int iBlock=nB*MAXN;
	double curSS=0;
	double deltaSS=0;
	double tSS=0;
	double t;
	double bestSS=0;
	int ni;
	int nj;
	double *deli=0;
	double *delj=0;


	*failure=true;

	ni=blocksizes[curBlock];

	curRowNo=B[IB(curBlock,xcur)];
	xri=X+curRowNo*k;
	xmi=blockMeans+curBlock*k;

	if (doWholeBlock) {
		deli=blockFactors+curBlock*k;
	}


	for (newBlk=0;newBlk<nB;newBlk++) {
		if (newBlk!=curBlock) {
			nj=blocksizes[newBlk];
			if (doWholeBlock)
				delj=blockFactors+newBlk*k;
			xmj=blockMeans+newBlk*k;

			tSS=blockSS[curBlock]+blockSS[newBlk];

			for (xnw=0;xnw<nj;xnw++) {
				newRowNo=B[IB(newBlk,xnw)];

				xrj=X+newRowNo*k;

				curSS=0;

				for (i=0;i<k;i++)
					vec[i]=(xrj[i]-xri[i])/ni; 
				if (doWholeBlock) {
					for (i=0;i<k;i++)
						vec[i]*=deli[i];

				}
				for (i=0;i<k;i++) {
					t=xmi[i]-gMean[i];
					t+=vec[i]; /* (mean-gmean)+(b-a)/ni */
					if (obScaled)
						curSS+=t*t/gSS[i];
					else 
						curSS+=t*t;
				}

				for (i=0;i<k;i++)
					vec[i]=(xri[i]-xrj[i])/nj;
				if (doWholeBlock) {
					for (i=0;i<k;i++)
						vec[i]*=delj[i];
				}
				for (i=0;i<k;i++) {
					t=xmj[i]-gMean[i];
					t+=vec[i];
					if (obScaled)
						curSS+=t*t/gSS[i];
					else 
						curSS+=t*t;
				}

				curSS=tSS-curSS;

				if (curSS>bestSS) {
					bestSS=curSS;
					*failure=false;
					delta=bestSS;
					*newBlock=newBlk;
					*xnew=xnw;
					return(delta);
				}
			}
		}

	}


	if (extraBlock) { /* See if extra rows can be swapped  */

		for (j=0;j<nEx;j++) {
			newRowNo=B[iBlock+j];

			xrj=X+newRowNo*k;

												 
			for (i=0;i<k;i++) {
				t=xrj[i]-xri[i];
				if (doWholeBlock)
					t*=deli[i];
				gMeanSwap[i]=gMean[i]+t/Nxb; /* corrected grand mean */
				vec[i]=t/ni;
			}

			curSS=0;				/* Get the difference between total SS's for gMean and gMeanSwap */
			for (i=0;i<nB;i++) {
				xmn=blockMeans+i*k;
				for (l=0;l<k;l++) {
					xrc[l]=xmn[l];
					if (i==curBlock)
						xrc[l]+=vec[l]; /* xrc is the current mean adjusted for the swap */
				}
				deltaSS=0;
				for (l=0;l<k;l++) {
					t=xrc[l]+gMean[l]-gMeanSwap[l];
					if (obScaled)
						deltaSS+=t*t/gSS[l];
					else 
						deltaSS+=t*t;	
				}
				curSS+=blockSS[i]-deltaSS;
			}


			if (curSS>bestSS) {
				bestSS=curSS;
				*failure=false;
				delta=curSS;
				*newBlock=nB;
				*xnew=j;
				return(delta);
			}
		}
	}


	return(delta);
}

/* findDeltaDpc ************************************************************************
| Calculates the delta for exchanging xold in the current block with xnew in a new block
*/ 

double findDeltaDpc(
	double *Tip,
	double *X,
	double *blockMeans,
	double *tX,
	double *tBlockMeans,
	double *vec,
	int		*B,
	int		nB,
	int		nEx,
	int		*blocksizes,
	double  *blockFactors,
	int		xcur,
	int		*xnew,
	int		curBlock,
	int		*newBlock,
	int		k,
	bool    *failure
	)
{
	double delta=0;
	double d;
	double d1;
	double d2;
	double Gi[3];
	double Gj[3];
	double M[3];
	double M1[3];
	double gi;
	double gj;
	double gij;
	double wii;
	double wij;
	double wji;
	double wjj;
	int    ni;
	int    nj;
	double *xri=tX; /* current point  */
	double *xrj=tX+k; /* new point -- to be made */
	double *xmi=tBlockMeans; /* current block mean  */
	double *xmj=tBlockMeans+k; /* new block mean -- to be made */
	double *deli;
	double *delj;
	int newBlk;
	int xnw;
	int K=(k*(k+1))/2;
	int i;
	int j;
	int l=0;
	int curRowNo;
	int newRowNo;
	int kc;
	int kn;
	int iBlock=nB*MAXN;

	ni=blocksizes[curBlock];

	kc=(ni<=k)?ni-1:k;

	curRowNo=B[IB(curBlock,xcur)];

	memcpy((void *)vec,X+curRowNo*k,k*sizeof(double));
	if (doWholeBlock) {
		deli=blockFactors+curBlock*k;
		for (i=0;i<kc;i++)
			vec[i]*=deli[i];
	}	
	transformVect(Tip+curBlock*K,vec,xri,kc); 

	memcpy((void *)vec,blockMeans+curBlock*k,k*sizeof(double));
	if (doWholeBlock) {
		deli=blockFactors+curBlock*k;
		for (i=0;i<kc;i++)
			vec[i]*=deli[i];
	}	
	transformVect(Tip+curBlock*K,vec,xmi,kc); 


	Gi[0]=(1+1/(double)ni);
	Gi[1]=1/(double)ni;
	Gi[2]=-(1-1/(double)ni);

	gi=0;
	for (l=0;l<kc;l++) {
		wii=xri[l]-xmi[l];
		gi+=wii*wii;
	}
	M[2]=gi;

	*failure=true;

	for (newBlk=0;newBlk<nB;newBlk++) {
		if (newBlk!=curBlock) {
			nj=blocksizes[newBlk];
			kn=(nj<=k)?nj-1:k;

			Gj[0]=(1+1/(double)nj);
			Gj[1]=1/(double)nj;
			Gj[2]=-(1-1/(double)nj);


			for (xnw=0;xnw<nj;xnw++) {
				newRowNo=B[IB(newBlk,xnw)];

				memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					deli=blockFactors+curBlock*k;
					for (i=0;i<kc;i++)
						vec[i]*=deli[i];
				}
				transformVect(Tip+curBlock*K,vec,xrj,kc); 

				gj=gij=0;
				for (l=0;l<kc;l++) {
					wii=xri[l]-xmi[l];
					wji=xrj[l]-xmi[l];
					gj+=wji*wji;
					gij+=wii*wii;
				}
				M[0]=gj;
				M[1]=gij;


				l=3;
				while(l--)
					M1[l]=Gi[l]+M[l];

				d1=-(M1[0]*M1[2]-M1[1]*M1[1]);
				d1=log(d1)/kc;

				memcpy((void *)vec,X+curRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					delj=blockFactors+newBlk*k;
					for (i=0;i<kn;i++)
						vec[i]*=delj[i];
				}
				transformVect(Tip+newBlk*K,vec,xri,kn);

				memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					delj=blockFactors+newBlk*k;
					for (i=0;i<kn;i++)
						vec[i]*=delj[i];
				}

				transformVect(Tip+newBlk*K,vec,xrj,kn);

				memcpy((void *)vec,blockMeans+newBlk*k,k*sizeof(double));
				if (doWholeBlock) {
					delj=blockFactors+newBlk*k;
					for (i=0;i<kn;i++)
						vec[i]*=delj[i];
				}

				transformVect(Tip+newBlk*K,vec,xmj,kn);


				gi=gj=gij=0;
				for (l=0;l<kn;l++) {
					wij=xri[l]-xmj[l];
					wjj=xrj[l]-xmj[l];
					gi+=wij*wij;
					gj+=wjj*wjj;
					gij+=wij*wjj;
				}
				M[0]=gi;
				M[1]=gij;
				M[2]=gj;

				l=3;
				while(l--)
					M1[l]=Gj[l]+M[l];

				d2=-(M1[0]*M1[2]-M1[1]*M1[1]);
				d2=log(d2)/kn;

				d=d1+d2;


				if ((d-delta)>deltaTol) {
					*failure=false;
					delta=d;
					*newBlock=newBlk;
					*xnew=xnw;
				}
			}
		}
	}

	if (extraBlock) { /* See if extra rows can be swapped  */
		for (j=0;j<nEx;j++) {
			newRowNo=B[iBlock+j];

			memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
			if (doWholeBlock) {
				deli=blockFactors+curBlock*k;
				for (i=0;i<kc;i++)
					vec[i]*=deli[i];
			}
			transformVect(Tip+curBlock*K,vec,xrj,kc); 

			gj=gij=0;
			for (l=0;l<kc;l++) {
				wii=xri[l]-xmi[l];
				wji=xrj[l]-xmi[l];
				gj+=wji*wji;
				gij+=wii*wii;
			}
			M[0]=gj;
			M[1]=gij;


			l=3;
			while(l--)
				M1[l]=Gi[l]+M[l];

			d=-(M1[0]*M1[2]-M1[1]*M1[1]);


			if ((d-delta)>deltaTol) {
				*failure=false;
				delta=d;
				*newBlock=nB;
				*xnew=j;
			}
		}
	}


	return(delta);
}

/* findDeltaDp ************************************************************************
| Calculates the delta for exchanging xold in the current block with xnew in a new block
*/ 

double findDeltaDp(
	double *Tip,
	double *X,
	double *tX,
	int		*B,
	int		nB,
	int     nEx,
	int		*blocksizes,
	double  *blockFactors,
	double  *vec,
	int		xcur,
	int		*xnew,
	int		curBlock,
	int		*newBlock,
	int		k,
	bool    *failure
	)
{
	double delta=0;
	double dx;
	double dy;
	double dyc;
	double dxy;
	double d;
	double d1;
	double d2;
	int    ni;
	int    nj;
	double *ypc=tX; /* transformed current point in current block -- to be replaced */
	double *xpc=tX+k; /* transformed new point in current block */
	double *yp=tX+2*k;  /* transformed current point in new block */
	double *xp=tX+3*k;	/* transformed new point in new block -- to be replaced */
	int newBlk;
	int xnw;
	int K=(k*(k+1))/2;
	int curRowNo;
	int newRowNo;
	int kc;
	int kn;
	double *deli=0;
	double *delj=0;
	int	   i;
	int	   j;
	int iBlock=nB*MAXN;

	ni=blocksizes[curBlock];

	kc=(ni>=k)?k:ni;



	curRowNo=B[IB(curBlock,xcur)];

	memcpy((void *)vec,X+curRowNo*k,k*sizeof(double));
	if (doWholeBlock) {
		deli=blockFactors+curBlock*k;
		for (i=0;i<kc;i++)
			vec[i]*=deli[i];
	}
	
	transformVect(Tip+curBlock*K,vec,ypc,kc); 

	dyc=0;
	for (i=0;i<kc;i++)
		dyc+=ypc[i]*ypc[i];
	
	*failure=true;

	for (newBlk=0;newBlk<nB;newBlk++) {
		if (newBlk!=curBlock) {
			nj=blocksizes[newBlk];

			for (xnw=0;xnw<nj;xnw++) {
				newRowNo=B[IB(newBlk,xnw)];

				memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					for (i=0;i<kc;i++)
						vec[i]*=deli[i];
				}
				transformVect(Tip+curBlock*K,vec,xpc,kc); 

				dx=dxy=0;
				for (i=0;i<kc;i++) {
					dx+=xpc[i]*xpc[i];
					dxy+=ypc[i]*xpc[i];
				}
				
					/* 1+delta */
				d1=(1+dx)*(1-dyc)+dxy*dxy; /* Replacing y with x in current block */
				d1=log(d1)/kc;

				kn=(nj>=k)?k:nj;

				memcpy((void *)vec,X+curRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					delj=blockFactors+newBlk*k;
					for (i=0;i<kn;i++)
						vec[i]*=delj[i];
				}
				transformVect(Tip+newBlk*K,vec,yp,kn);

				memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
				if (doWholeBlock) {
					for (i=0;i<kn;i++)
						vec[i]*=delj[i];
				}
				transformVect(Tip+newBlk*K,vec,xp,kn);

				dx=dy=dxy=0;
				for (i=0;i<kn;i++) {
					dx+=xp[i]*xp[i];
					dy+=yp[i]*yp[i];
					dxy+=xp[i]*yp[i];
				}
					/* 1+delta */
				d2=(1-dx)*(1+dy)+dxy*dxy;  /* Replacing x with y in new block */
				d2=log(d2)/kn;

				d=d1+d2;


				if ((d-delta)>deltaTol) {
					*failure=false;
					delta=d;
					*newBlock=newBlk;
					*xnew=xnw;
				}
			}
		}
	}

	if (extraBlock) { /* See if extra rows can be swapped  */
		for (j=0;j<nEx;j++) {
			newRowNo=B[iBlock+j];

			memcpy((void *)vec,X+newRowNo*k,k*sizeof(double));
			if (doWholeBlock) {
				for (i=0;i<kc;i++)
					vec[i]*=deli[i];
			}
			transformVect(Tip+curBlock*K,vec,xpc,kc); 

			dx=dxy=0;
			for (i=0;i<kc;i++) {
				dx+=xpc[i]*xpc[i];
				dxy+=ypc[i]*xpc[i];
			}
				/* 1+delta */
			d=(1+dx)*(1-dyc)+dxy*dxy; /* Replacing y with x in current block */
			d=log(d)/kc;

			if ((d-delta)>deltaTol) {
				*failure=false;
				delta=d;
				*newBlock=nB;
				*xnew=j;
			}
		}

	}

	return(delta);
}




/* transform **********************************************************************
| transfroms X and blockMeans to tX = X*Ti and tBlockMeans = tBlockMeans * Ti, 
|	using Tip which containts Ti'
*/

void transform(
	double *Tip,
	double *X,
	double *tX,
	double *blockMeans,
	double *tBlockMeans,
	int		N,
	int		k,
	int		nB
	)
{
	int i;
	int j;
	int l;
	double *pX;
	double *pXl;
	double *ptX;
	double *ptXl;
	double *pTip;
	double *pb;
	double *pbl;
	double *ptb;
	double *ptbl;

	memset((void *)tBlockMeans,0,nB*k*sizeof(double));
	memset((void *)tX,0,N*k*sizeof(double));

	for (i=0;i<N;i++) {
		pX=X+i*k;
		ptX=tX+i*k;
		pTip=Tip;
		for (j=0;j<k;j++) {
			pXl=pX;
			ptXl=ptX+j;
			for (l=0;l<=j;l++) {
				(*ptXl)+=(*(pXl++))*(*(pTip++));
			}
		}
	}
	
	for (i=0;i<nB;i++) {
		pb=blockMeans+i*k;
		ptb=tBlockMeans+i*k;
		pTip=Tip;
		for (j=0;j<k;j++) {
			pbl=pb;
			ptbl=ptb+j;
			for (l=0;l<=j;l++) {
				(*ptbl)+=(*(pbl++))*(*(pTip++));
			}
		}
		
	}

}





/* transformVect ****************************************************************
| transforms vec to tvec with Tip
*/

void transformVect(
	double *Tip,
	double *vec,
	double *tvec,
	int		k
	)
{
	int i;
	int j;
	double *pTip=Tip;
	double *pvec;
	double *ptvec;

	memset((void *)tvec,0,k*sizeof(double));
	ptvec=tvec;

	for (i=0;i<k;i++) {
		j=i+1;
		pvec=vec;
		while(j--) 
			(*ptvec)+=(*(pvec++))*(*(pTip++));	
		ptvec++;
	}


}


int ProgAllocate(
	int		**B,			/* nB x max(blocksizes) array of row numbers from X */
	double  **blockMeans,	/* nB x k matrix of block means */
	double  **tBlockMeans,  /* nB x k matrix of transformed block means */
	int	    **BlockArray,   /* N array of row numbers */
	double  **tX,			/* Transformed X */
	double	**T,			/* X'X=T'T, with T upper triangualar (has scale values on diagonal) */
	double  **Tip,			/* T inverse (multiplied by scale values) */
	double  **W,			/* k*(k+1)/2 scratch */
	double  **vec,			/* scratch 2*k element vector */
	double  **Sc,			/* scratch 2*k element vector */
	int		**rows,			/* scratch N element vector */
	int		N,				/* Number of rows in X */
	int		k,				/* Number of terms */
	int     Nxb,				/* Sum of block sizes */
	int		nB,				/* Number of blocks */
	bool    criterion,			/* If true, enlarge T and Tip */
	int		*blocksizes				/* Number of trials in the nB blocks */
)


{
	int K=(k*(k+1))/2;  
	int i;
	int	Nt=maxm(N,5*k);  /* Safety net for allocation in findDeltaBlockWhole() */
	int tBlock;

	MAXN=0;
	for (i=0;i<nB;i++)
		MAXN=maxm(MAXN,blocksizes[i]);
	nColumns=k;

	tBlock=nB*MAXN;
	if (extraBlock)
		tBlock+=N-Nxb;

	*B = (int *)malloc(tBlock * sizeof(int));
    if (!*B) return 4;
    
    *blockMeans = (double *)malloc(nB * k * sizeof(double));
    if (!*blockMeans) return 5;
    
    *tBlockMeans = (double *)malloc(nB * k * sizeof(double));
    if (!*tBlockMeans) return 5;
    
    *BlockArray = (int *)malloc(Nxb * sizeof(int));
    if (!*BlockArray) return 5;
    
    *tX = (double *)malloc(Nt * k * sizeof(double));
    if (!*tX) return 6;
    
    if (criterion>0) {
        *T = (double *)malloc(nB * K * sizeof(double));
        if (!*T) return 7;
        *Tip = (double *)malloc(nB * K * sizeof(double));
        if (!*Tip) return 8;
    }
    else {
        *T = (double *)malloc(K * sizeof(double));
        if (!*T) return 7;
        *Tip = (double *)malloc(K * sizeof(double));
        if (!*Tip) return 8;
    }
    
    *W = (double *)malloc(maxm(k*k,5*k) * sizeof(double)); /* 5*k needed in findDeltaBlockWhole() */
    if (!*W) return 9;
    
    *vec = (double *)malloc(2 * k * sizeof(double));
    if (!*vec) return 10;
    
    *Sc = (double *)malloc(2 * k * sizeof(double));
    if (!*Sc) return 11;
    
    *rows = (int *)malloc(maxm(N,Nxb) * sizeof(int));
    if (!*rows) return 12;


	return 0;
}

#ifdef NOTUSED
/* ProgDeallocate *********************************************************************
|	Frees space allocated by ProgAllocate
*/

void ProgDeallocate(
	int		*B,				/* nB x MAXN array */
	double  *blockMeans,	/* nB x k matrix of block means */
	double  *tBlockMeans,   /* nB x k matrix of transformed block means */
	double  *tX,			/* Transformed X plus block means */
	double	*T,			/* X'X=T'T, with T upper triangualar (has scale values on diagonal) */
	double  *Tip,			/* T inverse (multiplied by scale values) */
	double  *W,
	double  *vec,			/* scratch 2*k element vector */
	double  *Sc,			/* scratch 2*k element vector */
	int		*rows			/* scratch N element vector */
)
{
	if (B)
		R_Free(B);
	if (blockMeans)
		R_Free(blockMeans);
	if (tBlockMeans)
		R_Free(tBlockMeans);
	if (tX)
		R_Free(tX);
	if (T)
		R_Free(T);
	if (Tip)
		R_Free(Tip);
	if (W)
		R_Free(W);
	if (vec)
		R_Free(vec);
	if (Sc)
		R_Free(Sc);
	if (rows)
		R_Free(rows);

}

#endif

/* NoDupPermute *********************************************************************
|  Finds a permutation that will not cause duplicates in current block
*/

void NoDupPermuteB(
	int *rows,
	int N,
	int *B,
	int n,
	int bs
	)
{
	bool nodup=true;
	int i;
	int j;
	int curVal;

	repeat
		nodup=true;
		PermuteB(rows,N);
		for (i=0;i<n;i++) {
			curVal=B[i];
			for (j=0;j<bs-n;j++) {
				if (rows[j]==curVal) {
					nodup=false;
					break;
				}
			}
			if (!nodup)
				break;
		}
	until(nodup);

}


/* initializeB ***********************************************************************
|   Randomly assigns rows of X to blocks in B, recyciling if necessary.
*/

void initializeB(
	int *B,
	int *rows,
	int *irows,
	int N,
	int Nxb,
	int nB,
	int *blocksizes,
	bool firstRepeat

)
{
	int i;
	int j;
	int l;
	int bs;
	int iBlock=nB*MAXN;
	int t;
	int Nt=(initRows)?Nxb:N;

	// debug
	for (i=0;i<Nt;i++)
		rows[i]=i;

	if (initRows) { /* make irows the head of rows */
		for (i=0;i<Nxb;i++) {
			t=rows[i];
			rows[i]=irows[i];
			rows[irows[i]]=t;
		}
		if (!firstRepeat)	/* The user may have imposed a block strucure, so try it once */
			PermuteB(rows,Nxb);
	}
	else {
		PermuteB(rows,Nt);
	}



	for (i=0;i<nB*MAXN;i++)
		B[i]=-1;



	l=0;
	for (i=0;i<nB;i++) {
		bs=blocksizes[i];
		for (j=0;j<bs;j++) {
			if (l>=Nt) {
				l=0;
				NoDupPermuteB(rows,N,B+IB(i,0),j,bs);
			}
			B[IB(i,j)]=rows[l++];
		}
	}

	if (extraBlock) { /* Put the leftover rows in the extra block */
		for (i=l;i<Nt;i++)
			B[iBlock++]=rows[i];
	}


}

/* FillB ************************************************************************************
|  Fills B from BlockArray
*/

void FillB(
	int nB,
	int *B,
	int *blocksizes,
	int *BlockArray
)
{
	int i;
	int l=0;
	int bs;
	int j;

	for (i=0;i<nB;i++) {
		bs=blocksizes[i];
		for (j=0;j<bs;j++) {
			B[IB(i,j)]=BlockArray[l++]-1;
		}
	}


}

/* initializeBlockArray ***********************************************************************
|   Initialize BlockArray in case all designs are singular.
*/

void initializeBlockArray(
	int *rows,
	int *irows,
	int N,
	int Nxb,
	int nB,
	int *blocksizes,
	int *BlockArray

)
{
	int i;
	int j;
	int l;
	int m;
	int bs;
	int Nt=(initRows)?Nxb:N;

	if (initRows) {
		for (i=0;i<Nxb;i++)
			rows[i]=irows[i];
	}
	else {
		for (i=0;i<N;i++)
			rows[i]=i;
	}

	l=0;
	m=0;
	for (i=0;i<nB;i++) {
		bs=blocksizes[i];
		for (j=0;j<bs;j++) {
			if (l>=Nt) {
				l=0;
			}
			BlockArray[m++]=rows[l++]+1; /* In case there is no improvement when optimizing */
		}
	}

}


/* MeanAndSS **********************************************************************
|   updates mean and SS
*/

void MeanAndSS(
	double *x,
	double *mean,
	double *SS,
	int n,
	int k
)
{
	int np=n+1;
	int i;
	double delta;
	double meanInc;

	for (i=0;i<k;i++) {
		delta=x[i]-mean[i];
		meanInc=delta/np;
		if (1<n)
			SS[i]+=n*delta*meanInc;
		mean[i]+=meanInc;
	}
}

/* formBlockMeansOB *******************************************************************
|	Uses the row numbers in B to find block means.
|	Returns grand mean in gMean, the sum of squares of each block mean vector minus gMean, 
|   and the function value. 
|	NOTE: blockMeans and gMean are not scaled for whole block factors. This is done when
|		   they are used elsewhere. blockSS, however, are calculated from scaled means.
*/

double  formBlockMeansOB(
	double *blockMeans,
	double *X,
	int	   *B,
	int		k,
	int		nB,
	int     Nxb,
	int		*blocksizes,
	double  *blockFactors,
	double  *gMean,
	double  *gSS,
	double  *tolerance,
	double  *blockSS
)

{
	int i;
	int j;
	int l;
	int rowNo;
	double *xmi;
	double rn;
	double *pX;
	int  *pB;
	double fValue=0;
	double *deli=0;
	int    bs;
	int    n=0;
	double tTol=0;

	memset((void *)blockMeans,0,nB*k*sizeof(double));
	memset((void *)gMean,0,k*sizeof(double));
	memset((void *)gSS,0,k*sizeof(double));
	memset((void *)blockSS,0,nB*sizeof(double));



	for (i=0;i<nB;i++) {
		pB=B+i*MAXN;				/* pointer to first row no in block i */
		xmi=blockMeans+i*k;			/* pointer to first mean in block i */
		bs=blocksizes[i];
		for (j=0;j<bs;j++) {
			rowNo=*(pB++);
			pX=X+rowNo*k;			/* pointer to first element of row rowNo */
			MeanAndSS(pX,gMean,gSS,n++,k);
			for (l=0;l<k;l++) {
				xmi[l]+=*(pX++);
			}
		}

		if (doWholeBlock) {
			deli=blockFactors+i*k;
		}

		rn=(double)bs;
		for (l=0;l<k;l++) {
			if (doWholeBlock)
				xmi[l]*=deli[l]; /* the block means and gMean are scaled */
			xmi[l]/=rn;
		}

	}
	tTol=0;
	for (i=0;i<k;i++) {  /* form the grand SS */
		gSS[i]=gSS[i]/(n-1);
		tTol+=log(gSS[i]);
	}
	tTol/=k;
	*tolerance=exp(tTol); /* Used to decide upon zero in calling routine */

	for (i=0;i<nB;i++) {
		xmi=blockMeans+i*k;
		for (l=0;l<k;l++) {
			xmi[l]-=gMean[l];
			if (obScaled)
				blockSS[l]+=xmi[l]*xmi[l]/gSS[i];
			else 
				blockSS[l]+=xmi[l]*xmi[l];
		}
	}
	for (i=0;i<nB;i++) {
		fValue+=blockSS[i];
	}
	return(fValue);
}

	


/* formBlockMeans *******************************************************************
|	Uses the row numbers in B to find block means
*/

void formBlockMeans(
	double *blockMeans,
	double *X,
	int	   *B,
	int		k,
	int		nB,
	int		*blocksizes
)
{
	int i;
	int j;
	int l;
	int rowNo;
	double *xmi;
	double *pbmi;
	double rn;
	double *pX;
	int  *pB;

	memset((void *)blockMeans,0,nB*k*sizeof(double));


	for (i=0;i<nB;i++) {
		pB=B+i*MAXN;				/* pointer to first row no in block i */
		xmi=blockMeans+i*k;			/* pointer to first mean in block i */
		for (j=0;j<blocksizes[i];j++) {
			rowNo=*(pB++);
			pX=X+rowNo*k;			/* pointer to first element of row rowNo */
			pbmi=xmi;
			for (l=0;l<k;l++) {
				*(pbmi++)+=*(pX++);
			}
		}
		pbmi=xmi;
		rn=(double)blocksizes[i];
		for (l=0;l<k;l++) {
			*(pbmi++)/=rn;
		}
	}
}

	



/* BlockOptimize ********************************************************************
|	The optimizing function. 
|	Returns 13 in case the determinant is not positive
|
*/

#define designTol 1e-12


void BlockOptimize(
	double	*X, 
	int		nB,
	int		*blocksizes,
	double  *blockFactors,
	int		*B,
	double  *blockMeans,
	double  *tBlockMeans,
	double  *T,
	double  *Tip,
	double  *W,
	double  *tX,
	double  *vec,
	double  *Sc,
	int		*rows,
	int		*irows,
	int		N,
	int		Nxb,
	int		k,
	int		nEx,
	double  *D,
	double  *diagonality,
	int     *BlockArray,
	int		nRepeats,
	int     *iter,
	int     *error
)

{
	double	logDcrit=0;
	double  logDbest=0;
	double  delta;
	int		curBlock=0;
	int		newBlock=0;
	int		xcur;
	int		xnew;
	bool    singular;
	bool	exchanged;
	int		countSingular=0;
	int		nRepeatCounts=nRepeats;
	double  aVar=1;
	double  avVar=1;
	int	i;
	int j;


	initializeBlockArray(rows,irows,N,Nxb,nB,blocksizes,BlockArray);
	*iter=0;
	repeat{
		initializeB(B,rows,irows,N,Nxb,nB,blocksizes,nRepeatCounts==nRepeats);

		formBlockMeans(blockMeans,X,B,k,nB,blocksizes);

		logDcrit=reduceXtoT(X,T,B,blockMeans,k,nB,blocksizes,blockFactors,vec,Sc,&singular);


		if (!singular) {
			makeTiFromTB(Tip,T,W,&aVar,k);
			if (doWholeBlock) {
				repeat {
					exchanged=false;
					curBlock=0;
					repeat {
						for (xcur=0;xcur<blocksizes[curBlock];xcur++) {
							delta=findDeltaBlockWhole(X,Tip,W,blockMeans,B,nB,nEx,blocksizes,blockFactors,xcur,&xnew,
								curBlock,&newBlock,k);
							if (10>delta && delta>designTol) {
								exchangeBlockWhole(T,X,vec,blockMeans,B,blocksizes,blockFactors,xcur,xnew,curBlock,newBlock,nB,k);
								logDcrit+=log(1+delta);
								exchanged=true;
								makeTiFromTB(Tip,T,W,&aVar,k);
							}
							//R_CheckUserInterrupt();
						}
					} until(nB<=++curBlock);
				} until(!exchanged);

			}
			else {
				transform(Tip,X,tX,blockMeans,tBlockMeans,N,k,nB);
				repeat {
					exchanged=false;
					curBlock=0;
					repeat {
						for (xcur=0;xcur<blocksizes[curBlock];xcur++) {
							delta=findDeltaBlock(tX,tBlockMeans,B,nB,nEx,blocksizes,xcur,&xnew,curBlock,&newBlock,k);
							  /* poor starting designs cause numerical problems resulting in */
							  /* very large deltas */
							if (10>delta && delta>designTol) {
								/* one can insert formBlockMeans() and reduceXtoT() here to
								deal with numerical problems, but this seems not to be
								needed when large deltas are excluded, as above. */
								exchangeBlock(T,X,vec,blockMeans,B,blocksizes,xcur,xnew,curBlock,newBlock,nB,k);
								logDcrit+=log(1+delta);
								exchanged=true;
								makeTiFromTB(Tip,T,W,&aVar,k);
								transform(Tip,X,tX,blockMeans,tBlockMeans,N,k,nB);
							}
							//R_CheckUserInterrupt();
						}
					} until(nB<=++curBlock);
				} until(!exchanged);
			}

			if (logDcrit>logDbest) {
				(*iter)++;
				logDbest=logDcrit;
				avVar=aVar;
				j=0;
				for (i=0;i<nB*MAXN;i++) {
					if (B[i]!=-1) {
						BlockArray[j++]=B[i]+1;
					}
				}
			}
	

		}
		else 
			countSingular++;

		
	}until(!(--nRepeatCounts)); 

	if (countSingular==nRepeats)
		*error=13;
	else {

		*error=0;

		if (logDbest==0) {/*no improvement over original allocation */
			*error=22;
			*D=exp(logDcrit/(double)k)/(double)Nxb;
			*diagonality=0;
		}
		else {
			*D=exp(logDbest/(double)k)/(double)Nxb;
			*diagonality=1/((*D)*avVar*(double)Nxb);
		}
	}
}

/* BlockOptimizeOB ********************************************************************
|	Finds orthogonal blocks. 
|	Returns 13 in case the determinant is not positive
|
*/
#define SStol 1.0e-10

void BlockOptimizeOB(
	double	*X, 
	int		nB,
	int		*blocksizes,
	double  *blockFactors,
	int		*B,
	double  *blockMeans,
	double  *T,
	double  *W,
	double  *vec,
	double  *Sc,
	int		*rows,
	int		*irows,
	int		N,
	int		Nxb,
	int		k,
	int		nEx,
	double  *D,
	double  *diagonality,
	int     *BlockArray,
	int		nRepeats,
	int     *iter,
	int     *error
)


{
	double  delta;
	int		curBlock=0;
	int		newBlock=0;
	int		xcur=0;
	int		xnew=0;
	bool    singular;
	bool	exchanged;
	int		nRepeatCounts=nRepeats;
	int	i;
	int j;
	double *gMean=W;
	double *gSS=Sc;
	double *blockSS=W+k;
	double tolerance=0;
	double totalSS=0;
	double curSS=0;
	double bestSS=1e10;
	bool   failure;
	double logDet;

	initializeBlockArray(rows,irows,N,Nxb,nB,blocksizes,BlockArray);
	*iter=0;
	repeat{
		initializeB(B,rows,irows,N,Nxb,nB,blocksizes,nRepeatCounts==nRepeats);


		totalSS=formBlockMeansOB(blockMeans,X,B,k,nB,Nxb,blocksizes,blockFactors,gMean,gSS,&tolerance,blockSS);
		

		curSS=totalSS;
		tolerance*=SStol;
		repeat {
			exchanged=false;
			curBlock=0;
			repeat {
				for (xcur=0;xcur<blocksizes[curBlock];xcur++) {
					delta=findDeltaOB(X,blockMeans,vec,blockSS,B,nB,nEx,blocksizes,blockFactors,gMean,gSS,xcur,&xnew,
						curBlock,&newBlock,k,Nxb,&failure);

					if (!failure && delta>tolerance) {
						exchangeOB(X,vec,blockMeans,gMean,gSS,blockSS,B,blocksizes,blockFactors,xcur,xnew,
							curBlock,newBlock,nB,k,Nxb);
						curSS-=delta;
						exchanged=true;
					}
					//R_CheckUserInterrupt();
				}
			} until(nB<=++curBlock);
		} until(!exchanged);

		if (curSS<bestSS) {
			(*iter)++;
			bestSS=curSS;
			j=0;
			for (i=0;i<nB*MAXN;i++) {
				if (B[i]!=-1) {
					BlockArray[j++]=B[i]+1;
				}
			}
		}

	}until(!(--nRepeatCounts)); 

	FillB(nB,B,blocksizes,BlockArray);
	formBlockMeans(blockMeans,X,B,k,nB,blocksizes);
	logDet=reduceXtoT(X,T,B,blockMeans,k,nB,blocksizes,blockFactors,vec,Sc,&singular);
	*D=exp(logDet/(double)k)/(double)Nxb;


	*error=0;

	if (bestSS<tolerance)
		bestSS=0;
	*diagonality=bestSS;  /* A kluge */
}



/* BlockOptimizeDpc ********************************************************************
|	Optimizes Dpc 
|	Returns 13 in case some determinant is not positive, and 27 in case the design is singulare
|
*/



void BlockOptimizeDpc(
	double	*X, 
	int		nB,
	int		*blocksizes,
	double  *blockFactors,
	int		*B,
	double  *blockMeans,
	double  *tBlockMeans,
	double  *T,   /*This is a stack of triangular matrices. Fixed amount of memory is allocated for k(k+1)/2 each
				  event though when Transposed it use the triangular matrix may be smaller. */
	double  *Tip, /* is this */
	double  *W,
	double  *tX,
	double  *vec,
	double  *Sc,
	int     *rows,
	int     *irows,
	int		N,
	int		Nxb,
	int		nEx,
	int		k,
	double  *D,
	double  *Dp,
	int     *BlockArray,
	int		nRepeats,
	int		*iter,
	int     *error
)

{
	double	logDcrit=-1000;
	double  logDbest=-1000;
	double  logDet;
	double  delta;
	int		curBlock=0;
	int		newBlock=0;
	int		xcur;
	int		xnew;
	int		countSingular=0;
	bool    singular;
	bool    failure;
	bool	exchanged;
	int		nRepeatCounts=nRepeats;
	int		iterT=0;

	int	i;
	int j;
	int bs;

	initializeBlockArray(rows,irows,N,Nxb,nB,blocksizes,BlockArray);

	*iter=0;
	repeat{
		initializeB(B,rows,irows,N,Nxb,nB,blocksizes,nRepeatCounts==nRepeats);
			/*the determinants are scaled by the number of rows that are reduced, either blocksize or k */
		formBlockMeans(blockMeans,X,B,k,nB,blocksizes);
			/*returns an average determinant of blocks */
		logDcrit=reduceXtoTDpc(X,T,B,blockMeans,N,k,nB,blocksizes,blockFactors,vec,Sc,&singular);
		iterT=0;
		if (!singular) {
			makeTiFromTDpc(Tip,T,W,blocksizes,nB,-1,-1,k); /*Takes care of blocksizes internally */
			repeat {
				exchanged=false;
				curBlock=0;
				repeat {
					bs=blocksizes[curBlock];
					for (xcur=0;xcur<bs;xcur++) {
						delta=findDeltaDpc(Tip,X,blockMeans,tX,tBlockMeans,vec,B,nB,nEx,blocksizes,
							blockFactors,xcur,&xnew,curBlock,&newBlock,k,&failure);

						if (!failure && 1>delta && delta>designTol) {
							logDcrit+=delta;
							exchangeDpc(T,X,vec,blockMeans,B,blocksizes,blockFactors,xcur,xnew,curBlock,newBlock,nB,k);
							makeTiFromTDpc(Tip,T,W,blocksizes,nB,curBlock,newBlock,k);
						}
						//();
					}
				} until(nB<=++curBlock);
			} until(!exchanged || 100<iterT++); /*iter just in case */

			if (logDcrit>logDbest) {
				logDbest=logDcrit;
				(*iter)++;
				j=0;
				for (i=0;i<nB*MAXN;i++) {
					if (B[i]!=-1) {
						BlockArray[j++]=B[i]+1;
					}
				}
			}


		}
		else 
			countSingular++;

	}until(!(--nRepeatCounts)); 
	if (countSingular==nRepeats)
		*error=13;
	else {

		*error=0;
			/*Report D */
		FillB(nB,B,blocksizes,BlockArray);
		formBlockMeans(blockMeans,X,B,k,nB,blocksizes);
			/*This is used to get determinant of whole design, not average determinant of blocks */
		logDet=reduceXtoT(X,T,B,blockMeans,k,nB,blocksizes,blockFactors,vec,Sc,&singular);
		*D=exp(logDet/(double)k)/(double)Nxb;

		if (logDbest==-1000) {/*no improvement over original allocation */
			*error=22;
			*Dp=exp(logDcrit/nB);
		}
		else {
			*Dp=exp(logDbest/nB); /*(Prod[(|Mi/ni|)^(1/ki)])^(1/nB) */
		}
		if (singular)
			*error=27;

	}

}

/* BlockOptimizeDp ********************************************************************
|	Optimizes Dp 
|	Returns 13 in case some determinant is not positive, and 27 if the design is singular
|
*/




void BlockOptimizeDp(
	double	*X, 
	int		nB,
	int		*blocksizes,
	double  *blockFactors,
	int		*B,
	double  *blockMeans,
	double  *T,   /*This is a stack of triangular matrices. Fixed amount of memory is allocated for k(k+1)/2 each
				  event though when Transposed it use the triangular matrix may be smaller.*/
	double  *Tip, /*so is this */
	double  *W,
	double  *tX,
	double  *vec,
	double  *Sc,
	int     *rows,
	int     *irows,
	int		N,
	int		Nxb,
	int		nEx,
	int		k,
	double  *D,
	double  *Dp,
	int     *BlockArray,
	int		nRepeats,
	int		*iter,
	int     *error
)

{
	double	logDcrit=-1000;
	double  logDbest=-1000;
	double  logDet;
	double  delta;
	int		curBlock=0;
	int		newBlock=0;
	int		xcur;
	int		xnew;
	int		countSingular=0;
	bool    singular;
	bool    failure;
	bool	exchanged;
	int		nRepeatCounts=nRepeats;
	int		iterT=0;

	int	i;
	int j;
	int bs;

	initializeBlockArray(rows,irows,N,Nxb,nB,blocksizes,BlockArray);

	*iter=0;
	repeat{
		initializeB(B,rows,irows,N,Nxb,nB,blocksizes,nRepeatCounts==nRepeats);
			/*the determinants are scaled by the number of rows that are reduced, either blocksize or k */
		logDcrit=reduceXtoTDp(X,T,B,N,k,nB,blocksizes,blockFactors,vec,Sc,&singular);
		iterT=0;
		if (!singular) {
			makeTiFromTDp(Tip,T,W,blocksizes,nB,-1,-1,k); /*Takes care of blocksizes internally */
			repeat {
				exchanged=false;
				curBlock=0;
				repeat {
					bs=blocksizes[curBlock];
					for (xcur=0;xcur<bs;xcur++) {
						delta=findDeltaDp(Tip,X,tX,B,nB,nEx,blocksizes,blockFactors,vec,xcur,&xnew,curBlock,
							&newBlock,k,&failure);
						if (!failure && 1>delta && delta>designTol) {
							logDcrit+=delta;
							exchangeDp(T,X,vec,B,blocksizes,blockFactors,xcur,xnew,curBlock,newBlock,nB,k);
							makeTiFromTDp(Tip,T,W,blocksizes,nB,curBlock,newBlock,k);
						}
						//R_CheckUserInterrupt();
					}
				} until(nB<=++curBlock);
			} until(!exchanged || 100<iterT++); /*iter just in case */

			if (logDcrit>logDbest) {
				logDbest=logDcrit;
				(*iter)++;
				j=0;
				for (i=0;i<nB*MAXN;i++) {
					if (B[i]!=-1) {
						BlockArray[j++]=B[i]+1;
					}
				}
			}


		}
		else {
			countSingular++;
		}
	}until(!(--nRepeatCounts)); 
	if (countSingular==nRepeats)
		*error=13;
	else {

		*error=0;
			/*Report D */
		FillB(nB,B,blocksizes,BlockArray);
		formBlockMeans(blockMeans,X,B,k,nB,blocksizes);
		logDet=reduceXtoT(X,T,B,blockMeans,k,nB,blocksizes,blockFactors,vec,Sc,&singular);
		*D=exp(logDet/(double)k)/(double)Nxb;

		if (logDbest==-1000) {/*no improvement over original allocation */
			*error=22;
			*Dp=exp(logDcrit/nB);
		}
		else {
			*Dp=exp(logDbest/nB); /*(Prod[(|Mi/ni|)^(1/ki)])^(1/nB) */
		}
		if (singular)
			*error=27;

	}

}



void transposeMatrix(
	double *X,
	int N,  /* number of rows in transposed matrix */
	int k	/* number of columns in transposed matrix */
)	
{
    int     size = N*k-2;
	int		i=1;
	int		row;
	int		column;
	int		indx;
	double  temp;

	for(i=1;i<size;i++){
		indx = i;
		repeat
			column = indx/k; /* indx (column+k*row) of cell in transposed matrix */
			row = indx%k;
			indx = N*row +  column; /* the value to be swapped is in cell at indx in current matrix */
		solongas(indx < i); /* indx<i is a previously exchanged cell, */
							/* track its original contents and then exchange when found */

		if (indx >i) {
			temp=X[i];
			X[i] = X[indx];
			X[indx] = temp;
		}
	}
}

/* BlockOpt *****************************************************************************
| Command-line version of block optimization function
|
| Input: 
|   X - Matrix of model expanded design (N x k), without constant column
|   N - Number of rows
|   k - Number of terms (columns in X)
|   nB - Number of blocks
|   blocksizes - Array of block sizes
|   blockFactors - Matrix of block interaction factors (nB x k), or NULL
|   nRepeats - Number of optimization repeats
|   criterion - Optimization criterion (1-4)
|   initRows - Whether to use initial row assignments
|   irows - Initial row assignments (if initRows is true)
|
| Output:
|   Returns 0 on success, error code otherwise
|   Writes results to specified output files
*/
int BlockOpt(
    double *X,               // Input matrix (N x k)
    int N,                   // Number of rows
    int k,                   // Number of columns
    int nB,                 // Number of blocks
    int *blocksizes,        // Array of block sizes
    double *blockFactors,   // Block factors matrix (or NULL)
    int nRepeats,          // Number of repeats
    int criterion,         // Optimization criterion
    bool initRows,         // Whether to use initial rows
    int *irows            // Initial row assignments (or NULL)
) {
    double  D;
    double  Dp = 0.0;
    double  diagonality = 0.0;

    int     *BlockArray;
    int     *B;             // Block assignments
    double  *blockMeans;    // Block means
    double  *tBlockMeans;   // Transformed block means
    double  *tX;           // Transformed X
    double  *T;            // Upper triangular matrix
    double  *Tip;          // Inverse of T
    double  *W;            // Scratch matrix
    double  *vec;          // Scratch vector
    double  *Sc;           // Scratch vector
    int     *rows;         // Scratch array
    int     Nxb = 0;       // Sum of block sizes
    int     nEx = 0;       // Extra block size
    int    error;
    int     iter;

    // Calculate total block size
    for (int i = 0; i < nB; i++) {
        Nxb += blocksizes[i];
    }

    // Check if extra block needed
    extraBlock = false;
    if (Nxb < N) {
        extraBlock = true;
        nEx = N - Nxb;
    }

    // Set global flags
    doWholeBlock = (blockFactors != NULL);
    
    // Allocate memory
    error = ProgAllocate(&B, &blockMeans, &tBlockMeans, &BlockArray, &tX, &T, &Tip, 
                        &W, &vec, &Sc, &rows, N, k, Nxb, nB, criterion, blocksizes);
    
    if (error) {
        return error;
    }

    // Transpose input matrix for C row-major order
    transposeMatrix(X, N, k);
    if (doWholeBlock) {
        transposeMatrix(blockFactors, nB, k);
    }

    // Run optimization based on criterion
    obScaled = (criterion == 4);
    if (criterion == 3 || criterion == 4) {
        BlockOptimizeOB(X, nB, blocksizes, blockFactors, B, blockMeans, T, W, vec,
            Sc, rows, irows, N, Nxb, k, nEx, &D, &diagonality, BlockArray, nRepeats, &iter, &error);
    }
    else if (criterion == 1) {
        BlockOptimizeDpc(X, nB, blocksizes, blockFactors, B, blockMeans, tBlockMeans, T, Tip, W, tX, vec,
            Sc, rows, irows, N, Nxb, nEx, k, &D, &Dp, BlockArray, nRepeats, &iter, &error);
    }
    else if (criterion == 2) {
        BlockOptimizeDp(X, nB, blocksizes, blockFactors, B, blockMeans, T, Tip, W, tX, vec,
            Sc, rows, irows, N, Nxb, nEx, k, &D, &Dp, BlockArray, nRepeats, &iter, &error);
    }
    else {
        BlockOptimize(X, nB, blocksizes, blockFactors, B, blockMeans, tBlockMeans, T, Tip, W, tX, vec,
            Sc, rows, irows, N, Nxb, k, nEx, &D, &diagonality, BlockArray, nRepeats, &iter, &error);
    }

    // Write results to output files or print to stdout
    printf("Results:\n");
    printf("D criterion: %f\n", D);
    printf("Dp criterion: %f\n", Dp);
    printf("Diagonality: %f\n", diagonality);
    printf("Iterations: %d\n", iter);
    printf("Error code: %d\n", error);
    
    printf("\nBlock assignments:\n");
    for (int i = 0; i < Nxb; i++) {
        printf("%d ", BlockArray[i]);
    }
    printf("\n");

    // Clean up
    // Note: Memory allocated with R_alloc is automatically freed
    
    return error;
}






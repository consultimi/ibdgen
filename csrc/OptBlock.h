#define Imat(ROW,COL)  (COL)+(ROW)*nColumns-((ROW)*((ROW)+1))/2 /* Upper Triangular */
#define Iblock(ROW,COL) (COL)+(ROW)*nColumns					/* Index into block means */
#define IB(ROW,COL) (COL)+(ROW)*MAXN							/* Index into block list */

#define NOINVERSE 0
#define DOINVERSE 1
#define INVERSEONLY 2

int		MAXN=0;
int		nColumns=0;
bool    doWholeBlock=false; /* create block within block interactions */
bool    extraBlock=false;  /* True when candidate list is in an extra block */
bool	initRows=false;    /* True when initial design rows are specified */
bool    obScaled=false;    /* When true orthogonal blocks are scaled */

extern void transposeMatrix(double *X,int N,int k);

int BlockOpt(double *X,int N,int k,int nB,int *blocksizes,double *blockFactors,int nRepeats,int criterion,bool initRows,int *irows);

void BacksolveB(double *matrixXY,int	nTerms,int nColumns,int doInverse);
void Difference(double *vec,double *first,double *second,int	k);
void FillB(int nB,int *B,int *blocksizes,int *BlockArray);
void formBlockMeans(double *blockMeans,double *X,int *B,int k,int	nB,int *blocksizes);
void getRangeB(double *pMx,double *pMn,double *pX,int k);
void initializeB(int *B,int *rows,int *irows,int N,int Nxb,int nB,int *blocksizes,bool firstRepeat);
void initializeBlockArray(int *rows,int *irows,int N,int Nxb,int nB,int *blocksizes,int *BlockArray);
void MakeTi(double *T,double *Tip,int k);
void NoDupPermuteB(int *rows,int N,int *B,int n,int bs);
void PermuteB(int *a,int	n);
void RotateB(double	*vec,double	*tVec,double *matrixXY,int nTerms,int nColumns,
			double weight);
void transform(double *Tip,double *X,double *tX,double *blockMeans,double *tBlockMeans,
	int	N,int k,int	nB);
void transformVect(double *Tip,double *vec,double *tvec,int	k);
double transW(double *Tip,double *tVec,double *W,int k); 

int ProgAllocate(int **B,double **blockMeans,double **tBlockMeans,int **BlockArray,double **tX,double **T,		
	double **Tip,double **W,double **vec,double **Sc,int **rows,int N,int k,int Nxb,int nB,
	bool criterion,int *blocksize);

/* used only when Calloc() is used. */
/* void ProgDeallocate(int *B,double *blockMeans,double *tBlockMeans,double *tX,double *T, */
/*	double *Tip,double *W,double *vec,double *Sc,int *rows); */


/* The following are for the D criterion */
void BlockOptimize(double *X,int nB,int *blocksizes,double *blockFactors,int *B,double *blockMeans,
	double *tBlockMeans,double *T,double *Tip,double *W,double *tX,double *vec,
	double *Sc,int *rows, int *irows, int N,int Nxb,int k,int nEx,double *D,double *diagonality,
	int *BlockArray,int nRepeats,int *iter,int *error);
double reduceXtoT(double *X,double *T,int *B,double *blockMeans,int	k,int nB,
	int *blocksizes,double *blockFactors,double *vec,double *Sc,bool *singular);
void makeTiFromTB(double *Tip,double *T,double *W,double *avVar, int k);
double findDeltaBlock(double *tX,double *tBlockMeans,int	*B,int nB,int nEx,int *blocksizes,int xcur,
	int	*xnew,int curBlock,int *newBlock,int k);
double findDeltaBlockWhole(double *X,double *Tip,double *W,double *blockMeans,int *B,int nB,int nEx,int *blocksizes,
	double *blockFactors, int xcur,	int	*xnew,int curBlock,int *newBlock,int k);
void exchangeBlock(double *T,double *X,double *vec,double *blockMeans,
	int *B,int	*blocksizes,int	xcur,int xnew,int curBlock,int newBlock,int	nB,int k);
void exchangeBlockWhole(double *T,double *X,double *vec,double *blockMeans,int *B,int *blocksizes,double *blockFactors,
	int	xcur,int xnew,int curBlock,int newBlock,int	nB,int k);

/* The following are for the OB criterion */
void BlockOptimizeOB(double	*X,int	nB,int	*blocksizes,double  *blockFactors,int	*B,double  *blockMeans,
	double  *T,double  *W,double  *vec,double  *Sc,int *rows,
	int		*irows,int	N,int	Nxb,int	k,int	nEx,double  *D,double  *diagonality,int   *BlockArray,
	int		nRepeats,int  *iter,int  *error);
double  formBlockMeansOB(double *blockMeans,double *X,int   *B,int	k,int	nB,int    Nxb,int	*blocksizes,
	double  *blockFactors,double  *gMean,double *gSS,double *tolerance,double  *blockSS);
double findDeltaOB(double *X,double *blockMeans,double *vec,double *blockSS,int	*B,int	nB,int	nEx,
	int		*blocksizes,double  *blockFactors,double *gMean,double *gSS,int	xcur,int *xnew,int	curBlock,int	*newBlock,
	int		k,int	Nxb,bool   *failure);
void exchangeOB(double  *X,double  *vec,double  *blockMeans,double  *gMean,double *gSS,double  *blockSS,
	int		*B,int	*blocksizes,double *blockFactors,int	xcur,int xnew,int	curBlock,
	int		newBlock,int nB,int	k,int	Nxb);
void MeanAndSS(double *x,double *mean,double *SS,int n,int k);



/* The following are for the Dpc criterion */
void BlockOptimizeDpc(double *X,int nB,int *blocksizes,double *blockFactors,int *B,double *blockMeans,
	double *tBlockMeans,double *T,double *Tip,double *W,double *tX,double *vec,double *Sc,int *rows,int *irows,
	int N,int Nxb,int nEx,int k,double *D,double *Dp,int *BlockArray,int nRepeats,int *iter,int *error);
double reduceXtoTDpc(double *X,double *T,int *B,double *blockMeans,int N,int k,int nB,
	int *blocksizes,double *blockFactors,double *vec,double *Sc,bool *singular);
void makeTiFromTDpc(double *Tip,double *T,double	*W,int *blocksizes,int nB,int curBlock,int newBlock,int k);
double findDeltaDpc(double *Tip,double *X,double *blockMeans,double *tX,double *tBlockMeans,
	double *vec,int	*B,int nB,int nEx,int *blocksizes,double *blockFactors,int xcur,int *xnew,int curBlock,
	int	*newBlock,int k,bool *failure);
void exchangeDpc(double *T,double *X,double *vec,double *blockMeans,int *B,int *blocksizes,double *blockFactors,
	int	xcur,int xnew,int curBlock,int newBlock,int	nB,int k);

/* The following duplicate the above for the Dp criterion */
void BlockOptimizeDp(double *X,int nB,int *blocksizes,double *blockFactors,int *B,double *blockMeans,
	double *T,double *Tip,double *W,double *tX,double *vec,double *Sc,int *rows,int *irows,
	int N,int Nxb,int nEx,int k,double *D,double *Dp,int *BlockArray,int nRepeats,int *iter,int *error);
double reduceXtoTDp(double *X,double *T,int *B,int N,int k,int nB,
	int *blocksizes,double *blockFactors,double *vec,double *Sc,bool *singular);
void makeTiFromTDp(double *Tip,double *T,double	*W,int *blocksizes,int nB,int curBlock,int newBlock,int k);
double findDeltaDp(double *Tip,double *X,double *tX,int	*B,int nB,int nEx,int *blocksizes,double *blockFactors,
		double *vec,int xcur,int *xnew,int curBlock,int *newBlock,int	k,bool *failure);
void exchangeDp(double *T,double *X,double *vec,int *B,int *blocksizes,double *blockFactors,
	int	xcur,int xnew,int curBlock,int newBlock,int	nB,int k);


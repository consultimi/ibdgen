use nalgebra::{DMatrix, DVector};


/*/* initializeBlockArray ***********************************************************************
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

} */
fn initialize_block_array(rows: &mut Vec<u8>, irows: Vec<u8>, n: u8, n_xb: u8, n_b: u8, block_sizes: Vec<u8>, block_array: &mut Vec<u8>, init_rows: bool) -> Result<(), String> {

    if init_rows {
        for i in 0..n_xb {
            rows[i as usize] = irows[i as usize];
        }
    } else {
        for i in 0..n {
            rows[i as usize] = i;
        }
    }

    let mut l = 0;
    let mut m = 0;
    let n_t = if init_rows { n_xb } else { n };
    for i in 0..n_b {
        let bs = block_sizes[i as usize];
        for _ in 0..bs {
            if l >= n_t {
                l = 0;
            }
            block_array[m as usize] = rows[l as usize] + 1;
            m += 1;
            l += 1;
        }
    }
    Ok(())
}

/* PermuteB **********************************************************
|	Randomly pemutes the n integers in a[] using the Fike
|	algorithm.  See Fike, "A permutation generation method"  The Computer
|	Journal, 18-1, Feb 75, 21-22.
*/
fn permute_b(a: &mut Vec<u8>, n: u8) -> Result<(), String> {

    let rnd = rand::random::<f64>();
    for i in 1..n {
        let j = (((1 + i) as f64) * rnd) as i32;
        let temp = a[j as usize];
        a[j as usize] = a[i as usize];
        a[i as usize] = temp;
    }
    
    Ok(())
}

fn no_dup_permute_b(b: &mut DMatrix<f64>, rows: &mut Vec<u8>, n: u8, i: u8, bs: u8) -> Result<(), String> {
/*bool nodup=true;
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
	until(nodup); */
    let mut nodup = true;
    loop {
        nodup = true;
        permute_b(rows, n)?;
        for i in 0..n {
            let cur_val = b[(i * n + i) as usize];
            for j in 0..bs-n {
                if rows[j as usize] as f64 == cur_val {
                    nodup = false;
                    break;
                }
            }
            if !nodup {
                break;
            }
        }
        if nodup {
            break;
        }
    }

    Ok(())
}

fn initialize_b(b: &mut DMatrix<f64>, rows: &mut Vec<u8>, irows: Vec<u8>, n: u8, n_xb: u8, n_b: u8, 
    block_sizes: Vec<u8>, n_repeat_counts: u8, init_rows: bool, first_repeat: bool, max_n: u8, extra_block: bool) -> Result<(), String> {

    let mut i = 0;
    let mut j = 0;
    let mut l = 0;
    let mut bs = 0;
    let mut i_block = n_b * max_n;
    let mut t = 0;
    let mut n_t = if init_rows { n_xb } else { n };

    for i in 0..n_t {
        rows[i as usize] = i;
    }

    if init_rows {
        for i in 0..n_xb {
            t = rows[i as usize];
            rows[i as usize] = irows[i as usize];
            rows[irows[i as usize] as usize] = t;
        }
        if !first_repeat {
            permute_b(rows, n_xb)?;
        }
    } else {
        permute_b(rows, n_t)?;
    }

    /*	for (i=0;i<nB*MAXN;i++)
		B[i]=-1; */
    for i in 0..n_b * max_n {
        b[i as usize] = -1.0;
    }

    /*l=0;
	for (i=0;i<nB;i++) {
		bs=blocksizes[i];
		for (j=0;j<bs;j++) {
			if (l>=Nt) {
				l=0;
				NoDupPermuteB(rows,N,B+IB(i,0),j,bs);
			}
			B[IB(i,j)]=rows[l++];
		}
	} */
    let mut l = 0;
    for i in 0..n_b {
        bs = block_sizes[i as usize];
        for j in 0..bs {
            if l >= n_t {
                l = 0;
                no_dup_permute_b(rows, n, b, i, bs)?;
            }
            b[(i * max_n + j) as usize] = rows[l as usize] as f64;
            l += 1;
        }
    }

    /*	if (extraBlock) { /* Put the leftover rows in the extra block */
		for (i=l;i<Nt;i++)
			B[iBlock++]=rows[i];
	} */
    if extra_block {
        for i in l..n_t {
            b[(i_block + i) as usize] = rows[i as usize] as f64;
        }
    }

    Ok(())

}
// optimize determinant over all blocks using d-criterion
fn block_optimize(x: DMatrix<f64>, n_b: u8, block_sizes: Vec<u8>, block_factors: Option<DMatrix<f64>>, n_repeats: u8, 
        init_rows: bool, rows: &mut Vec<u8>, irows: Vec<u8>, n: u8, n_xb: u8, k: u8) -> Result<(), String> {
    /*
	int		*B, //  an nB x MAXN array of rows of X allocated to the nB blocks, -1 is used to fill blocks up to MAXN. N-nB*MAXN ints are added to hold the extra block when N is greater then the sum of the block sizes.
	double  *blockMeans, // nB x k matrix of block means
	double  *tBlockMeans, // transformed block means
	double  *T, // X'X=T'T, with T upper triangualar (has scale values on diagonal)
	double  *Tip, // T inverse (multiplied by scale values)
	double  *W, // k*(k+1)/2 scratch matrix
	double  *vec, // scratch 2*k element vector
	double  *Sc, // scratch 2*k element vector
	int		*rows, // rows to use
	int		*irows,
	int		N,
	int		Nxb,
	int		k,
	int		nEx,
	double  *D,
	double  *diagonality, 
	int     *BlockArray, // array of row numbers from X arranged in nB blocks
	int     *iter,
	int     *error
     */

    let mut b: DMatrix<f64> = DMatrix::zeros(n_b as usize, *block_sizes.iter().max().unwrap() as usize);
    let mut block_means: DMatrix<f64> = DMatrix::zeros(n_b as usize, k as usize);
    let mut t_block_means: DMatrix<f64> = DMatrix::zeros(n_b as usize, k as usize);
    let mut t: DMatrix<f64> = DMatrix::zeros(n as usize, n as usize);
    let mut tip: DMatrix<f64> = DMatrix::zeros(n as usize, k as usize);
    let mut w: DMatrix<f64> = DMatrix::zeros(k as usize * (k as usize + 1) / 2, 1);
    let mut vec: DVector<f64> = DVector::zeros(2 * k as usize);
    let mut sc: DVector<f64> = DVector::zeros(2 * k as usize);
    let mut block_array: Vec<u8> = vec![0; n_b as usize * (*block_sizes.iter().max().unwrap()) as usize];

    // b is a matrix of block factors. ncols is max(blocksizes)
    let mut b: DMatrix<f64> = DMatrix::zeros(n_b as usize, *block_sizes.iter().max().unwrap() as usize);
    initialize_block_array(rows, irows, n, n_xb, n_b, block_sizes, &mut block_array, init_rows)?;

    for _ in 0..n_repeats {

    }

    Ok(())
}

pub fn opt_block(x_i: DMatrix<f64>, init_rows: bool, rows: Vec<u8>, n_b: u8, block_sizes: Vec<u8>, do_whole_block: bool, 
    block_factors: Option<DMatrix<f64>>, n_repeats: u8, criterion: u8) -> Result<(), String> {
    
    let mut x = x_i.clone();
    let n = x.nrows() as u8;
    let k = x.ncols() as u8;

    let mut n_xb: u8 = 0;
    for i in 0..n_b {
        n_xb += block_sizes[i as usize];
    }

    let mut extra_block = false;
    if n_xb < n {
        extra_block = true;
    }

    x = x.transpose();

    if do_whole_block == true {
        if let Some(block_factors) = block_factors {
            let mut block_factors = block_factors.transpose();
        }
    }

    Ok(())
}


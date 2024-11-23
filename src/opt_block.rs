use nalgebra::{DMatrix, DVector, dmatrix};




fn reduce_x_to_t(
    block_data: &mut BlockData,
    vec: &mut DVector<f64>,
    sc: &mut DVector<f64>,
    do_whole_block: bool,
) -> (f64, bool) {
    let mut log_det = 0.0;
    let mut singular = false;

    //dbg!(&block_data.x, &block_data.x.clone().qr().r());
    let (p_mx, p_mn) = sc.as_mut_slice().split_at_mut(block_data.k as usize);
    p_mx.fill(f64::NEG_INFINITY);
    p_mn.fill(f64::INFINITY);

    block_data.t.fill(0.0);

    let mut tvec = DVector::zeros(block_data.k as usize);

    for i in 0..block_data.n_b {
        dbg!(&i, &block_data.max_n, &block_data.n_b);
        let p_b = block_data.b.view((i as usize * block_data.max_n as usize, 0), (block_data.max_n as usize, 1));
        let x_mi = block_data.block_means.row(i as usize);
        let deli = if do_whole_block {
            Some(block_data.block_factors.as_ref().unwrap().view((i as usize * block_data.k as usize, 0), (block_data.k as usize, block_data.k as usize)))
        } else {
            None
        };

        for j in 0..block_data.block_sizes[i as usize] {
            let row_no = p_b[(j as usize, 0)];
            let x_ri = block_data.x.row(row_no as usize);

            vec.iter_mut().zip(x_ri.iter().zip(x_mi.iter()))
                .for_each(|(v, (&x, &m))| *v = x - m);

            if let Some(ref deli) = deli {
                vec.component_mul_assign(deli);
            }

            get_range_b(p_mx, p_mn, vec.as_slice(), block_data.k as usize);
            rotate_b(vec, &mut tvec, &mut block_data.t, 1.0);
        }
    }

    for i in 0..block_data.k as usize {
        let r = (p_mx[i] + p_mn[i]) / 2.0;
        let t_val = block_data.t[(i, i)];
        if t_val <= 0.0 || t_val < r * 1e-10 {
            singular = true;
            return (0.0, singular);
        }
        log_det += t_val.ln();
    }

    (log_det, singular)
}

fn get_range_b(p_mx: &mut [f64], p_mn: &mut [f64], vec: &[f64], k: usize) {
    for i in 0..k {
        p_mx[i] = p_mx[i].max(vec[i]);
        p_mn[i] = p_mn[i].min(vec[i]);
    }
}

fn rotate_b(vec: &DVector<f64>, tvec: &mut DVector<f64>, matrix_xy: &mut DMatrix<f64>, weight: f64) {
    // Implementation of RotateB function goes here
    // This function would need to be implemented using nalgebra operations
    // It's a complex function and its implementation would depend on the specific requirements of your project
}

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
fn initialize_block_array(block_data: &mut BlockData, block_array: &mut Vec<u8>) -> Result<(), String> {

    if block_data.init_rows {
        for i in 0..block_data.n_xb {
            block_data.rows[i as usize] = block_data.irows[i as usize];
        }
    } else {
        for i in 0..block_data.n {
            block_data.rows[i as usize] = i;
        }
    }

    let mut l = 0;
    let mut m = 0;
    let n_t = if block_data.init_rows { block_data.n_xb } else { block_data.n };
    for i in 0..block_data.n_b {
        let bs = block_data.block_sizes[i as usize];
        for _ in 0..bs {
            if l >= n_t {
                l = 0;
            }
            block_array[m as usize] = block_data.rows[l as usize] + 1;
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
fn permute_b(a: &mut DVector<u8>, n: u8) -> Result<(), String> {

    let rnd = rand::random::<f64>();
    for i in 1..n {
        let j = (((1 + i) as f64) * rnd) as i32;
        let temp = a[j as usize];
        a[j as usize] = a[i as usize];
        a[i as usize] = temp;
    }
    
    Ok(())
}

fn no_dup_permute_b(block_data: &mut BlockData, little_n: u8, bs: u8) -> Result<(), String> {
    loop {
        dbg!(&block_data.rows);
        let mut nodup = true;
        permute_b(&mut block_data.rows, block_data.n)?;
        for i in 0..little_n {
            dbg!("{:?} {:?}", bs, little_n);
            let cur_val = block_data.b[(i * block_data.n_t + i) as usize];
            for j in 0..(bs - little_n) {
                if block_data.rows[j as usize] as f64 == cur_val {
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

fn form_block_means(block_data: &mut BlockData) {
    println!("block_data.b: {:?}", block_data.b);
    // divide block_data.b into block_data.n_b equal sized blocks of max_n rows
    // block_data.b is a n_b x max_n matrix of row indices from block_data.x
    // block_data.block_means is a n_b x k matrix of block means

    //block_data.block_means.fill(0.0);

    let mut i = 0;
    block_data.block_means.row_iter_mut().for_each(|mut block| {
        dbg!(&block);
        let block_row = block_data.b.row(i as usize);
        let out: Vec<_> = block_row.iter().map(|&row_index| {
            block_data.x.row(row_index as usize)
        }).collect();
        //dbg!(&out_mat.row_sum());
        
        dbg!(&DMatrix::from_rows(&out).row_sum());
        block.copy_from(&(DMatrix::from_rows(&out).row_sum() / block_data.block_sizes[i as usize] as f64));
        i += 1;
    });
    dbg!(&block_data.block_means);

}

fn initialize_b(block_data: &mut BlockData, first_repeat: bool) -> Result<(), String> {

    //let mut i = 0;
    let mut j = 0;
    let mut l = 0;
    let mut bs = 0;
    let mut i_block = block_data.n_b * block_data.max_n;
    let mut t = 0;


    for i in 0..block_data.n_t {
        block_data.rows[i as usize] = i;
    }

    if block_data.init_rows {
        for i in 0..block_data.n_xb {
            t = block_data.rows[i as usize];
            block_data.rows[i as usize] = block_data.irows[i as usize];
            block_data.rows[block_data.irows[i as usize] as usize] = t;
        }
        if !first_repeat {
            permute_b(&mut block_data.rows, block_data.n_xb)?;
        }
    } else {
        permute_b(&mut block_data.rows, block_data.n_t)?;
    }

    /*	for (i=0;i<nB*MAXN;i++)
		B[i]=-1; */
    for i in 0..block_data.n_b * block_data.max_n {
        block_data.b[i as usize] = -1.0;
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
    for i in 0..block_data.n_b {
        bs = block_data.block_sizes[i as usize];
        for j in 0..bs {
            if l >= block_data.n_t {
                l = 0;
                no_dup_permute_b(block_data, j, bs)?;
            }
            block_data.b[(i * block_data.max_n + j) as usize] = block_data.rows[l as usize] as f64;
            l += 1;
        }
    }

    /*	if (extraBlock) { /* Put the leftover rows in the extra block */
		for (i=l;i<Nt;i++)
			B[iBlock++]=rows[i];
	} */
    if block_data.extra_block {
        for i in l..block_data.n_t {
            block_data.b[(i_block + i) as usize] = block_data.rows[i as usize] as f64;
        }
    }

    Ok(())

}

// optimize determinant over all blocks using d-criterion
fn block_optimize(block_data: &mut BlockData, n_repeats: u8) -> Result<(), String> {
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

    let mut tip: DMatrix<f64> = DMatrix::zeros(block_data.n as usize, block_data.k as usize);
    let mut w: DMatrix<f64> = DMatrix::zeros(block_data.k as usize * (block_data.k as usize + 1) / 2, 1);
    let mut vec: DVector<f64> = DVector::zeros(2 * block_data.k as usize);
    let mut sc: DVector<f64> = DVector::zeros(2 * block_data.k as usize);
    let mut block_array: Vec<u8> = vec![0; block_data.n_b as usize * (*block_data.block_sizes.iter().max().unwrap()) as usize];

    // b is a matrix of block factors. ncols is max(blocksizes)
    initialize_block_array(block_data, &mut block_array)?;

    for _ in 0..n_repeats {
        initialize_b(block_data,  false)?;
        dbg!(&block_data.b);
        form_block_means(block_data);
        let (log_det, singular) = reduce_x_to_t(block_data, &mut vec, &mut sc, false);
        dbg!(&log_det, &singular);
    }

    Ok(())
}

#[derive(Debug)]
struct BlockData {
    x: DMatrix<f64>,    // x is the input matrix, typically it'll be a dummy coded design matrix (diag is 1, off-diag is 0, first row all 0)
    b: DMatrix<f64>,    // b is a matrix of block factors. ncols is max(blocksizes)
    block_means: DMatrix<f64>, // block_means is a vector of block means
    t_block_means: DMatrix<f64>, // t_block_means is a vector of transformed block means
    t: DMatrix<f64>,    // t is a matrix of transformed data. It's upper triangular and has scale values on the diagonal
    max_n: u8,          // max_n is the maximum block size
    n: u8,              // n is the number of rows in x
    k: u8,              // k is the number of columns in x
    n_t: u8,            // n_t is the number of rows to use. If init_rows is true, n_t = n_xb, otherwise n_t = n
    extra_block: bool,   // extra_block is true if there are extra rows in x that are not used in the blocks
    n_xb: u8,           // n_xb is the number of rows in x that are used in the blocks
    n_b: u8,            // n_b is the number of blocks
    block_sizes: Vec<u8>, // block_sizes is a vector of block sizes
    n_repeat_counts: u8, // n_repeat_counts is the number of repeats
    init_rows: bool,     // init_rows is true if the rows are initialized from irows
    rows: DVector<u8>,  // rows is a vector of row indices
    irows: DVector<u8>,
    block_factors: Option<DMatrix<f64>>
}

impl BlockData {
    fn new() -> Self {
        BlockData { 
            x: DMatrix::zeros(0, 0),
            b: DMatrix::zeros(0, 0), 
            block_means: DMatrix::zeros(0, 0), 
            t_block_means: DMatrix::zeros(0, 0), 
            t: DMatrix::zeros(0, 0),
            max_n: 0,
            n: 0,
            k: 0,
            n_t: 0,
            extra_block: false,
            n_xb: 0,
            n_b: 0,
            block_sizes: vec![0; 0],
            n_repeat_counts: 0,
            init_rows: false,
            rows: DVector::zeros(0),
            irows: DVector::zeros(0),
            block_factors: None }
    }

    fn rows_vec_length(&self) -> usize {
        std::cmp::max(self.n, self.n_xb) as usize
    }
}

pub fn opt_block(x_i: DMatrix<f64>, rows: Option<Vec<u8>>, n_b: u8, block_sizes: Vec<u8>, do_whole_block: bool, 
    block_factors: Option<DMatrix<f64>>, n_repeats: u8, _criterion: u8) -> Result<(), String> {
    
    let mut block_data = BlockData::new();
    block_data.x = x_i.clone();
    block_data.n = x_i.nrows() as u8;
    block_data.k = x_i.ncols() as u8;
    block_data.block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
    block_data.t_block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
    block_data.t = DMatrix::zeros(block_data.n as usize, block_data.n as usize);
    block_data.n_repeat_counts = n_repeats;
    block_data.max_n = *block_sizes.iter().max().unwrap();
    block_data.n_xb = block_data.block_sizes.iter().sum();
    block_data.n_b = n_b;
    block_data.block_sizes = block_sizes;
    if let Some(rows) = rows {
        block_data.init_rows = true;
        block_data.rows = rows.into();
    } else {
        block_data.init_rows = false;
        block_data.rows = DVector::zeros(block_data.rows_vec_length());
    }
    block_data.b = DMatrix::zeros(n_b as usize, block_data.max_n as usize);


    block_data.extra_block = if block_data.n_xb < block_data.n { true } else { false };
    //block_data.x = block_data.x.transpose();
    block_data.n_t = if block_data.init_rows { block_data.n_xb } else { block_data.n };
    
    /*
    if do_whole_block == true {
        if let Some(block_factors) = block_factors {
            block_data.block_factors = Some(block_factors.transpose());
        }
    } */

    //dbg!(&block_data);

    block_optimize(&mut block_data, n_repeats)?;

    Ok(())
}



mod tests {
    use super::*;

    fn dm7choose3() -> DMatrix<u8> {
        nalgebra::dmatrix![
            0,0,0,0,0,0,0,
            0,1,0,0,0,0,0,
            0,0,1,0,0,0,0,
            0,0,0,1,0,0,0,
            0,0,0,0,1,0,0,
            0,0,0,0,0,1,0,
            0,0,0,0,0,0,1
        ]
    }

    #[test]
    fn test_opt_block() {
        let x = dm7choose3();
        let result = opt_block(x.cast::<f64>(), None, 7, vec![3, 3,3,3,3,3,3], true, None, 1, 0);
        assert!(result.is_ok());
    }
}


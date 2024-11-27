use nalgebra::{DMatrix, DVector, dmatrix};
use pretty_print_nalgebra::*;


const DESIGN_TOL: f64 = 1.0e-10;

fn reduce_x_to_t(
    block_data: &mut BlockData,
    vec: &mut DVector<f64>,
    sc: &mut DVector<f64>,
    do_whole_block: bool,
) -> (f64, bool) {
    let mut log_det = 0.0;
    let mut p_mx: Vec<f64> = vec![-1e16; block_data.k as usize];
    let mut p_mn: Vec<f64> = vec![1e16; block_data.k as usize];

    let mut i = 0;
    let mut block_means = block_data.block_means.clone();
    let b_matrix = block_data.b.clone();
    for block in block_means.row_iter_mut() {
        //dbg!(&block);
        let block_row = b_matrix.row(i as usize);
        for &row_index in block_row.iter() {
            //println!("Calculating diff for row: {:#?} and block: {:#?}", block_data.x.row(row_index as usize), &block );
            //pretty_print!(&block_data.x.row(row_index as usize));
            println!("row : {}",pretty_print!(&block_data.x.row(row_index as usize)));
            println!("block means: {}",pretty_print!(&block));
            let diff = block_data.x.row(row_index as usize) - &block;
            println!("diff: {}",pretty_print!(&diff.transpose()));
            get_range_b(&mut p_mx, &mut p_mn, &diff.transpose(), block_data.k as usize);
            // do the rotation. 
            rotate_b(block_data, &diff.transpose(), 1.0);
            println!("t not pretty: {:?}", &block_data.t);
            
            
            //diff
        }
        //dbg!(&block);
        
        //dbg!(&DMatrix::from_rows(&out).row_sum());
        //block.copy_from(&(DMatrix::from_rows(&out).row_sum() / block_data.block_sizes[i as usize] as f64));
        i += 1;
        //println!("block_data.t: {:#?}", &block_data.t);
    };
    //let upper_tri = block_data.x.clone().lu().u();
    block_data.t = block_data.t.transpose();
    println!("t: {}", pretty_print!(&block_data.t));
    //println!("upper_tri: {:?}", &upper_tri);
    log_det = block_data.t.upper_triangle().determinant().log10();
    /*
    println!("log_det: {}", log_det);
    let mut t_cnt = 0;
    for i in 0..block_data.k {
        let r = (p_mx[i as usize] + p_mn[i as usize]) / 2.0;
        let t = block_data.t[t_cnt as usize];
        if t <= 0.0 || t < r * 1e-10 {
            return (0.0, true);
        }
        log_det += t.log10();
        t_cnt += (block_data.k - i) as usize;
    } */

    (log_det, false)
}


fn get_range_b(p_mx: &mut Vec<f64>, p_mn: &mut Vec<f64>, vec: &DVector<f64>, k: usize) {
    for i in 0..k {
        p_mx[i] = p_mx[i].max(vec[i]);
        p_mn[i] = p_mn[i].min(vec[i]);
    }
}

const TOLROT: f64 = 1.0e-12;

// this function does not map to the diagonal index of the upper triangular matrix
//fn calc_index(i: usize, nc: usize) -> usize {
//    i+i*nc-(i*(i+1))/2
//}

fn calc_index(i: usize, nc: usize) -> usize {
    i * (nc+1)
}

fn rotate_b(block_data: &mut BlockData, vec: &DVector<f64>, starting_weight: f64) {
    let mut skip: bool = false;
    let mut weight = starting_weight;
    // clone the diff vector
    let mut t_vec = vec.clone();
    //println!("vec: {:?}", &vec);
    //println!("t_vec: {:?}", &t_vec);
    let mut k_index = 0;
    for i in 0..block_data.k {
        if skip == true {
            break;
        }

        if t_vec[i as usize] == 0.0 {
            continue;
        }

        // d points to the corresponding index in the t (upper triangular) matrix
        k_index = calc_index(i as usize, block_data.k as usize);
        println!("i: {}, k_index: {}", i, k_index);
        let d = block_data.t[(k_index) as usize];
        let dp = d + weight * t_vec[i as usize] * t_vec[i as usize];
        if dp.abs() < TOLROT {
            continue;
        }
        //dbg!(&k_index);
        block_data.t[(k_index) as usize] = dp;

        let c = d / dp;
        let s = weight * t_vec[i as usize] / dp;

        if d == 0.0 {
            skip = true;
            weight = 0.0;
            continue;
        } else {
            weight *= c;
        }

        k_index += 1;
        for j in (i+1)..block_data.k {
            dbg!(&k_index);
            let r = block_data.t[k_index];
            block_data.t[k_index] = s * t_vec[j as usize] + c * r;
            t_vec[j as usize] -= t_vec[i as usize] * r;
            k_index += 1;
        }
    }
}

fn make_ti_from_tb(block_data: &mut BlockData) -> Result<f64, String> {
    println!("block_data.t: {}", pretty_print!(&block_data.t));
    block_data.t_inv = block_data.t.upper_triangle().clone().try_inverse().ok_or("Failed to invert upper triangular matrix")?;
    block_data.t_inv.fill_lower_triangle_with_upper_triangle();
    println!("ti: {}", pretty_print!(&block_data.t_inv));
    //let tip = ti.lower_triangle();
    
    // set the lower triangular part of block_data.t to tip
    //block_data.t.set_lower_triangle(tip);
    for i in 0..block_data.k {
        for j in 0..block_data.k {
            if j > i {
                block_data.t_inv[(i * block_data.k + j) as usize] = block_data.t_inv[(i * block_data.k + j) as usize];
            }
        }
    }

    
    println!("block_data.t: {}", pretty_print!(&block_data.t));
    println!("block_data.t_inv: {}", pretty_print!(&block_data.t_inv));

    Ok(())
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
        //dbg!(&block_data.rows);
        let mut nodup = true;
        permute_b(&mut block_data.rows, block_data.n)?;
        for i in 0..little_n {
            //dbg!("{:?} {:?}", bs, little_n);
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
    for mut block in block_data.block_means.row_iter_mut() {
        //dbg!(&block);
        let block_row = block_data.b.row(i as usize);
        let out: Vec<_> = block_row.iter().map(|&row_index| {
            block_data.x.row(row_index as usize)
        }).collect();
        //dbg!(&out_mat.row_sum());
        
        //dbg!(&DMatrix::from_rows(&out).row_sum());
        block.copy_from(&(DMatrix::from_rows(&out).row_sum() / block_data.block_sizes[i as usize] as f64));
        i += 1;
    };
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

fn find_delta_block(block_data: &mut BlockData, xcur: u8, cur_block: u8) -> Result<f64, String> {
    Ok(0.0)
}

fn exchange_block(block_data: &mut BlockData, xcur: u8, cur_block: u8, exchanged: &mut bool) -> Result<(), String> {
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
    let mut var = 0.0;
    // b is a matrix of block factors. ncols is max(blocksizes)
    initialize_block_array(block_data, &mut block_array)?;

    for repeat_num in 0..n_repeats {
        initialize_b(block_data,  repeat_num == 0)?;
        //dbg!(&block_data.b);
        form_block_means(block_data);
        let (log_det, singular) = reduce_x_to_t(block_data, &mut vec, &mut sc, false);
        if singular {
            return Err("Singular matrix".to_string());
        } else {
            make_ti_from_tb(block_data)?;
            /* transform **********************************************************************
| transfroms X and blockMeans to tX = X*Ti and tBlockMeans = tBlockMeans * Ti, 
            |	using Tip which containts Ti'
            */
            //transform(block_data)?;
            block_data.t_x = block_data.x.clone() * block_data.t_inv.clone();
            block_data.t_block_means = block_data.block_means.clone() * block_data.t_inv.clone();

            //println!("420 block_data.t_x: {}", pretty_print!(&block_data.t_x));
            //println!("421block_data.t_block_means: {}", pretty_print!(&block_data.t_block_means));
            loop {  
                let mut exchanged = false;
                let cur_block = 0;
                loop {
                    for xcur in 0..block_data.block_sizes[cur_block as usize] {
                        let delta = find_delta_block(block_data, xcur, cur_block)?;
                        if delta > 10.0 && delta > DESIGN_TOL {
                            exchange_block(block_data, xcur, cur_block, &mut exchanged)?;
                        }
                    }
                }
                dbg!(&log_det, &singular);
            }
        }
    }

    Ok(())
}

#[derive(Debug)]
struct BlockData {
    x: DMatrix<f64>,    // x is the input matrix, typically it'll be a dummy coded design matrix (diag is 1, off-diag is 0, first row all 0)
    t_x: DMatrix<f64>,  // x_t is the transpose of x
    b: DMatrix<f64>,    // b is a matrix of block factors. ncols is max(blocksizes)
    block_means: DMatrix<f64>, // block_means is a vector of block means
    t_block_means: DMatrix<f64>, // t_block_means is a vector of transformed block means
    t: DMatrix<f64>,    // t is a matrix of transformed data. It's upper triangular and has scale values on the diagonal
    t_inv: DMatrix<f64>, // t_inv is the inverse of t
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
            t_x: DMatrix::zeros(0, 0),
            b: DMatrix::zeros(0, 0), 
            block_means: DMatrix::zeros(0, 0), 
            t_block_means: DMatrix::zeros(0, 0), 
            t: DMatrix::zeros(0, 0),
            t_inv: DMatrix::zeros(0, 0),
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

pub fn opt_block(x_i: DMatrix<f64>, rows: Option<Vec<u8>>, n_b: u8, block_sizes: Vec<u8>,  n_repeats: u8) -> Result<(), String> {
    
    let mut block_data = BlockData::new();
    block_data.x = x_i.clone();
    block_data.n = x_i.nrows() as u8;
    block_data.k = x_i.ncols() as u8;
    block_data.block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
    block_data.t_block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
    block_data.t = DMatrix::zeros(block_data.k as usize, block_data.k as usize);
    block_data.t_inv = DMatrix::zeros(block_data.k as usize, block_data.k as usize);
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
    let test_matrix: DMatrix<u8> = nalgebra::dmatrix![
            1,2,3;
            4,5,6;
            7,8,9;
        ];
    println!("test_matrix: {}", pretty_print!(&test_matrix.upper_triangle()));
    println!("test_matrix: {:#?}", &test_matrix);

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
        let result = opt_block(x.cast::<f64>(), None, 7, vec![3, 3,3,3,3,3,3], 1);
        assert!(result.is_ok());
    }
}


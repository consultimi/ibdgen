use nalgebra::{DMatrix, DVector, SVector};
use pretty_print_nalgebra::*;
use anyhow::*;
use derive_builder::Builder;
const DESIGN_TOL: f64 = 1.0e-10;

const DEBUG: bool = false;

macro_rules! debug_println {
    ($($arg:tt)*) => {
        if DEBUG {
            println!($($arg)*);
        }
    };
}

#[derive(Builder, Debug)]
#[builder(build_fn(error = "anyhow::Error"))]
struct BlockData {

    /* x is the input matrix, typically it'll be a dummy coded design matrix (diag is 1, off-diag is 0, first row all 0) */
    x: DMatrix<f64>,   

    /* t_x is the transpose of x */
    t_x: DMatrix<f64>,  

    /* b is a matrix of block factors. ncols is max(blocksizes) */
    b: DMatrix<i32>,    

    /* block_means is a matrix of block means */
    block_means: DMatrix<f64>, 

    /* t_block_means is a matrix of transformed block means */
    t_block_means: DMatrix<f64>, 

    /* t is a matrix of transformed data. It's upper triangular and has scale values on the diagonal */
    t: DMatrix<f64>, 

    /* t_inv is the inverse of t */
    t_inv: DMatrix<f64>,

    /* max_n is the maximum block size */
    max_n: u8,

    /* n is the number of rows in x */
    n: u8,

    /* k is the number of columns in x */
    k: u8, 

    /* n_t is the number of rows to use. If init_rows is true, n_t = n_xb, otherwise n_t = n */
    n_t: u8,

    /* n_xb is the number of rows in x that are used in the blocks */
    n_xb: u8,           

    /* n_b is the number of blocks */
    n_b: u8,           
    
    #[builder(default = "RandomType::Uniform")]
    random_type: RandomType, // random_type is the type of randomization
    
    #[builder(setter(custom))]
    block_sizes: Vec<u8>, // block_sizes is a vector of block sizes

    /* rows is a vector of row indices */
    rows: DVector<u8>, 

    /* prohibited_pairs is a vector of prohibited pairs */
    #[builder(default = "vec![]")]
    prohibited_pairs: Vec<(u8, u8)>,

}

impl BlockData {


    fn reduce_x_to_t(&mut self) -> (f64, bool) {
        let mut p_mx: Vec<f64> = vec![-1e16; self.k as usize];
        let mut p_mn: Vec<f64> = vec![1e16; self.k as usize];
    
        // initialise T
        self.t.fill(0.0);
        let mut i = 0;
        let mut block_means = self.block_means.clone();
        let b_matrix = self.b.clone();
        for block in block_means.row_iter_mut() {
            //dbg!(&block);
            let block_row = b_matrix.row(i as usize);
            for &row_index in block_row.iter() {
                let diff = self.x.row(row_index as usize) - &block;
                get_range_b(&mut p_mx, &mut p_mn, &diff.transpose(), self.k as usize);
                self.rotate_b(&diff.transpose(), 1.0);
            }
            i += 1;
        };
        let log_det = self.t.determinant().ln();
        (log_det, false)
    }

    fn rotate_b(&mut self, vec: &DVector<f64>, starting_weight: f64) {

        debug_println!("rotate_b called with t: {}", pretty_print!(&self.t));
        let mut skip: bool = false;
        let mut weight = starting_weight;
        // clone the diff vector
        let mut t_vec = vec.clone();
        //debug_println!("vec: {:?}", &vec);
        //debug_println!("block_data.t: {}", pretty_print!(&block_data.t.transpose()));
        //let mut k_index = 0;
        for i in 0..self.k {
            if skip == true {
                break;
            }
    
            if t_vec[i as usize] == 0.0 {
                continue;
            }
    
            // d points to the corresponding index in the t (upper triangular) matrix
            let mut k_index = calc_index(i as usize, self.k as usize);
            //debug_println!("i: {}, k_index: {}", i, k_index);
            let d = self.t[(k_index) as usize];
            //debug_println!("d: {}", d);
            let dp = d + weight * t_vec[i as usize] * t_vec[i as usize];
            if dp.abs() < TOLROT {
                continue;
            }
            //dbg!(&k_index);
            self.t[(k_index) as usize] = dp;
            //debug_println!("t after set: {}", pretty_print!(&block_data.t));
            let c = d / dp;
            let s = weight * t_vec[i as usize] / dp;
    
            if d == 0.0 {
                skip = true;
                weight = 0.0;
                //continue;
            } else {
                weight *= c;
            }
    
            k_index += 1;
            for j in (i+1)..self.k {
                let r = self.t[k_index];
                self.t[k_index] = s * t_vec[j as usize] + c * r;
                t_vec[j as usize] -= t_vec[i as usize] * r;
                k_index += 1;
            }
        }
    }

    
    fn make_ti_from_tb(&mut self) -> Result<f64> {
        debug_println!("block_data.t (make_ti_from_tb): {}", pretty_print!(&self.t));
        //block_data.t_inv = block_data.t.upper_triangle().clone().try_inverse().ok_or(anyhow!("Failed to invert upper triangular matrix"))?;
        let mut t_inv = self.t.clone();
        t_inv.fill_diagonal(1.0);
        let mut t_inv = t_inv.try_inverse().ok_or(anyhow!("Failed to invert upper triangular matrix"))?;

        let diag = self.t.diagonal().apply_into(|x| *x = 1.0 / *x);
        t_inv.set_diagonal(&diag);
        
        let scale_vec = t_inv.diagonal().map(|x| x.sqrt());
        t_inv.fill_diagonal(1.0);

        
        for (i, mut row) in t_inv.row_iter_mut().enumerate() {
            row *= scale_vec[i];
        } 

        let rowtally: Vec<f64> = t_inv.column_iter().map(|x| x.map(|x| x.powi(2)).sum()).collect::<Vec<_>>();

        let mut a_var = rowtally.iter().map(|x| x.ln()).sum::<f64>() / self.k as f64;
        a_var = a_var.exp();
        //debug_println!("coltally (make_ti_from_tb): {:?}", &coltally);
        self.t_inv = t_inv;
        Ok(a_var)
    }

    fn initialize_block_array(&mut self, block_array: &mut Vec<u8>)  {

        for i in 0..self.n {
            self.rows[i as usize] = i;
        }
    
        let mut l = 0;
        let mut m = 0;
        let n_t = self.n;
        for i in 0..self.n_b {
            let bs = self.block_sizes[i as usize];
            for _ in 0..bs {
                if l >= n_t {
                    l = 0;
                }
                block_array[m as usize] = self.rows[l as usize] + 1;
                m += 1;
                l += 1;
            }
        }
    }

    fn no_dup_permute_b(&mut self, offset: u8, little_n: u8, bs: u8) -> Result<()> {
        //eprintln!("no_dup_permute_b called with little_n: {}, bs: {}, block_data.rows: {}, offset: {}", little_n, bs, &self.rows, offset);
        //debug_println!("b: {}", pretty_print!(&block_data.b));
        loop {
            let mut nodup = true;
            permute_b(&mut self.rows, self.n, self.random_type).map_err(|e| anyhow!("Failed to permute rows: {}", e))?;
            //eprintln!("rows after permute: {}", pretty_print!(&self.rows));
            for i in 0..little_n {
                //let index = offset * self.max_n + i;
                let cur_val = self.b[(offset as usize, i as usize)] as i32;
                for j in 0..(bs - little_n) {
                    if self.rows[j as usize] as i32 == cur_val {
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
    
    fn form_block_means(&mut self) {
        // divide block_data.b into block_data.n_b equal sized blocks of max_n rows
        // block_data.b is a n_b x max_n matrix of row indices from block_data.x
        // block_data.block_means is a n_b x k matrix of block means
        let mut i = 0;
        for mut block in self.block_means.row_iter_mut() {
            let block_row = self.b.row(i as usize);
            let out: Vec<_> = block_row.iter().map(|&row_index| {
                self.x.row(row_index as usize)
            }).collect();
            block.copy_from(&(DMatrix::from_rows(&out).row_sum() / self.block_sizes[i as usize] as f64));
            i += 1;
        };
        debug_println!("block_means inside form_block_means: {}", pretty_print!(&self.block_means));
    }

    
    fn initialize_b(&mut self) -> Result<()> {
        // 1. Initialize sequential numbers 0..n_t
        for i in 0..self.n_t {
            self.rows[i as usize] = i;
        }

        // 2. Randomly permute these numbers
        permute_b(&mut self.rows, self.n_t, self.random_type)?;

        // 3. Initialize b matrix with -1 (empty slots)
        for i in 0..self.n_b * self.max_n {
            self.b[((i / self.max_n) as usize, (i % self.max_n) as usize)] = -1;
        }

        let mut l = 0;  // Index into rows array
        
        // 4. For each block
        'block_loop: for i in 0..self.n_b {
            let bs = self.block_sizes[i as usize];  // Size of current block
            let mut retry_count = 0;
            
            // 5. Fill each position in the block
            let mut j = 0;
            while j < bs {
                // 6. If we've used all available points, reset and reshuffle
                if l >= self.n_t {
                    l = 0;
                    permute_b(&mut self.rows, self.n_t, self.random_type)?;
                }

                // 7. Get next candidate point
                let candidate = self.rows[l as usize];
                let mut is_valid = true;

                // 8. Check if candidate violates any prohibited pairs
                for k in 0..j {
                    let existing_point = self.b[(i as usize, k as usize)] as u8;
                    if self.prohibited_pairs.iter().any(|&(a, b)| 
                        (candidate == a && existing_point == b) || 
                        (candidate == b && existing_point == a)
                    ) {
                        is_valid = false;
                        break;
                    }
                }

                // 9a. If valid, place point and move to next position
                if is_valid {
                    self.b[(i as usize, j as usize)] = candidate as i32;
                    l += 1;
                    j += 1;
                    retry_count = 0;
                } 
                // 9b. If invalid, try next point
                else {
                    l += 1;
                    
                    const MAX_RETRIES: u8 = 10;
                    // 10. If we've tried all points, restart block
                    if l >= self.n_t {
                        retry_count += 1;
                        if retry_count >= MAX_RETRIES {
                            return Err(anyhow!("Unable to find valid configuration"));
                        }
                        
                        // Reset and try again
                        l = 0;
                        permute_b(&mut self.rows, self.n_t, self.random_type)?;
                        for k in 0..j {
                            self.b[(i as usize, k as usize)] = -1;
                        }
                        j = 0;
                    }
                }
            }
        }
        Ok(())
    }

    
    fn find_delta_block(&mut self, xcur: u8, xnew: &mut u8, cur_block: u8, new_block: &mut u8) -> Result<f64> {
        const DELTA_TOL: f64 = 1e-12;  // Minimum improvement threshold
        let mut delta = 0.0;  // Tracks best improvement found
        
        // Get current point's row number and block size
        let cur_row_no = self.b[(cur_block as usize, xcur as usize)] as usize;
        let ni = self.block_sizes[cur_block as usize];

        let fi = self.t_x.row(cur_row_no);
        let fmi = self.t_block_means.row(cur_block as usize);

        // Loop through all blocks except current
        for i in 0..self.n_b {
            if i != cur_block {
                let nj = self.block_sizes[i as usize];
                
                // Check if exchange would create prohibited combination in either block
                let would_violate_constraints = |candidate_row: u8| {
                    // Check target block - get all points except the one we're exchanging
                    let target_block_points: Vec<u8> = self.b.row(i as usize)
                        .iter()
                        .take(nj as usize)
                        .filter(|&&x| x >= 0)  // Filter out empty slots (-1)
                        .map(|&x| x as u8)
                        .filter(|&x| x != candidate_row)  // Exclude the point we're exchanging
                        .collect();

                    // Check current block - get all points except the one being moved
                    let current_block_points: Vec<u8> = self.b.row(cur_block as usize)
                        .iter()
                        .take(self.block_sizes[cur_block as usize] as usize)
                        .filter(|&&x| x >= 0)
                        .map(|&x| x as u8)
                        .filter(|&x| x != cur_row_no as u8)  // Exclude the point we're moving out
                        .collect();

                    // Check if current point would violate constraints with target block
                    for &point in &target_block_points {
                        if self.prohibited_pairs.iter().any(|&(a, b)| 
                            (cur_row_no as u8 == a && point == b) || 
                            (cur_row_no as u8 == b && point == a)
                        ) {
                            return true;
                        }
                    }

                    // Check if candidate point would violate constraints with current block
                    for &point in &current_block_points {
                        if self.prohibited_pairs.iter().any(|&(a, b)| 
                            (candidate_row == a && point == b) || 
                            (candidate_row == b && point == a)
                        ) {
                            return true;
                        }
                    }
                    
                    false
                };

                let g_vec: SVector<f64, 3> = SVector::from_vec(vec![
                    (ni + nj) as f64 / (ni * nj) as f64, 
                    1.0, 
                    0.0
                ]);

                let fmj = self.t_block_means.row(i as usize);
                let diff = fmj - fmi;
                let mut mi_vec: SVector<f64, 3> = SVector::from_vec(vec![
                    diff.component_mul(&diff).sum(), 
                    0.0, 
                    0.0
                ]);

                // Try exchanging with each point in candidate block
                for j in 0..nj {
                    //let row_no = b_transpose[(i * self.max_n + j) as usize] as usize;
                    let row_no = self.b[(i as usize, j as usize)] as usize;
                    // Skip if exchange would create prohibited combination
                    if self.prohibited_pairs.len() > 0 && would_violate_constraints(row_no as u8) {
                        //println!("proposed swap from block {} ({}) idx {} val {} with block {} ({}) idx {} val {} skipped", cur_block, self.b.row(cur_block as usize), xcur, cur_row_no, i, self.b.row(i as usize), j, row_no);
                        continue;
                    }
                    //println!("i is {}, j is {}", i, j);

                    let fj = self.t_x.row(row_no);
                    mi_vec[1] = (fmj - fmi).component_mul(&(fj - fi)).sum();
                    mi_vec[2] = (fj - fi).component_mul(&(fj - fi)).sum();


                    
                    // Combine geometric and moment terms
                    
                    
                    // Combine geometric and moment terms
                    let mlivec = g_vec + mi_vec;
                    let d = -(1.0 + mlivec[0] * mlivec[2] - mlivec[1] * mlivec[1]);

                    if (d - delta) > DELTA_TOL {
                        delta = d;
                        *new_block = i;
                        *xnew = j;
                    }
                }
            }
        }

        Ok(delta)
    }

    
    fn exchange_block(&mut self, xcur: u8, xnew: u8, cur_block: u8, new_block: &mut u8) -> Result<()> {
        let row_no_i = self.b[(cur_block as usize, xcur as usize)] as usize;
        let ni = self.block_sizes[cur_block as usize];

        let x_clone = self.x.clone();
        let xri = x_clone.row(row_no_i);
        let xmi = self.block_means.row(cur_block as usize);
        //debug_println!("xri: {}\nxmi: {}\nrowNoi: {}", pretty_print!(&xri), pretty_print!(&xmi), row_no_i);

        let row_no_j = self.b[(*new_block as usize, xnew as usize)] as usize;
        // Handle normal block exchange case

        let xrj = x_clone.row(row_no_j);
        let xmj = self.block_means.row(*new_block as usize);
        debug_println!("xmi: {}\nxmj: {}\nxri: {}\nxrj: {}\nrowNoj: {}", pretty_print!(&xmi), pretty_print!(&xmj), pretty_print!(&xri), pretty_print!(&xrj), row_no_j);
        let nj = self.block_sizes[*new_block as usize];
        let c = (ni + nj) as f64 / (ni * nj) as f64;

        // vec = xmj - xmi
        let mut vec = (xmj - xmi).transpose();
        self.rotate_b(&vec, 1.0);
        vec -= (xrj - xri).transpose();
        self.rotate_b(&vec, -1.0);
        vec = (xrj - xri).transpose();
        self.rotate_b(&vec, 1.0 - c);
        
        // Update block means
        for i in 0..self.k {
            //let idx = cmi_from_rmi((cur_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            //self.block_means[idx as usize] += newsum;
            self.block_means[(cur_block as usize, i as usize)] += (xrj[i as usize] - xri[i as usize]) / ni as f64;
            //let idx = cmi_from_rmi((*new_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            self.block_means[(*new_block as usize, i as usize)] += (xri[i as usize] - xrj[i as usize]) / nj as f64;
        }

        //println!("new_block: {}, xnew: {}, b before exchange: {}", *new_block, xnew, pretty_print!(&self.b));

        self.b[(*new_block as usize, xnew as usize)] = row_no_i as i32;
        self.b[(cur_block as usize, xcur as usize)] = row_no_j as i32;

        //println!("cur_block: {}, xcur: {}, b after exchange: {}", cur_block, xcur, pretty_print!(&self.b));

        Ok(())
    }

    fn block_optimize(&mut self, n_repeats: u8) -> Result<BlockResult> {

        let mut block_array: Vec<u8> = vec![0; self.n_b as usize * (*self.block_sizes.iter().max().unwrap()) as usize];
        let mut best_log_det = 0.0;
        let mut best_block_array = DMatrix::zeros(self.n_b as usize, self.max_n as usize);
        let mut xnew = 0;
        let mut new_block = 0;
        let mut av_var = 0.0;
        // b is a matrix of block factors. ncols is max(blocksizes)
        self.initialize_block_array(&mut block_array);
    
        for repeat_num in 0..n_repeats {
            println!("REPEAT NUMBER: {}", repeat_num);

            //rintln!("b before initialize_b: {}", pretty_print!(&self.b));
            self.initialize_b().map_err(|e| anyhow!("Failed to initialize b: {}", e))?;
            //println!("b after initialize_b: {}", pretty_print!(&self.b));

            self.form_block_means();
            let (mut log_det, singular) = self.reduce_x_to_t();
            if singular {
                return Err(anyhow!("Singular matrix"));
            } else {
                av_var = self.make_ti_from_tb().map_err(|e| anyhow!("Failed to make ti from tb: {}", e))?;
                self.transform();
                loop {  
                    let mut exchanged = false;
                    for cur_block in 0..self.n_b {
                        for xcur in 0..self.block_sizes[cur_block as usize] {
                            self.try_exchange_block(&mut xnew, &mut new_block, &mut av_var, &mut log_det, &mut exchanged, cur_block, xcur)?;
                        }
                    }
                    if !exchanged {
                        break;
                    }
    
                    if log_det > best_log_det {
                        best_log_det = log_det;
                        best_block_array = self.b.clone().try_cast::<usize>().unwrap();
                    }
                }
            }
        }
    
        debug_println!("best_log_det: {}", best_log_det);
        //println!("block_data.k: {}", self.k);
        //println!("block_data.n_xb: {}", self.n_xb);
        // 	*D=exp(logDbest/(double)k)/(double)Nxb;
        let best_d = (best_log_det / self.k as f64).exp() / self.n_xb as f64;
        let best_diagonality = 1.0 / (best_d * av_var * self.n_xb as f64);
        let best_coincidence = CoincidenceMatrix::from_block_array(&best_block_array);
        let best_block_array = BlockArray::from_block_array(&best_block_array);

        Ok(BlockResult { best_log_det, best_block_array, best_d, best_diagonality, best_coincidence })
    }

    fn try_exchange_block(&mut self, xnew: &mut u8, new_block: &mut u8, av_var: &mut f64, log_det: &mut f64, exchanged: &mut bool, cur_block: u8, xcur: u8) -> Result<(), Error> {
        debug_println!("BEING LOOP xcur: {}, curBlock: {}, newBlock: {}", xcur, cur_block, new_block);
        let delta = self.find_delta_block(xcur, xnew, cur_block, new_block).map_err(|e| anyhow!("Failed to find delta block: {}", e))?;
        debug_println!("delta: {}", delta);
        Ok(if delta < 10.0 && delta > DESIGN_TOL {
    
            debug_println!("t before exchange: {}", pretty_print!(&self.t));
            debug_println!("t_inv before exchange: {}", pretty_print!(&self.t_inv));
            self.exchange_block(xcur, *xnew, cur_block, new_block).map_err(|e| anyhow!("Failed to exchange block: {}", e))?;
            debug_println!("t_inv after exchange: {}", pretty_print!(&self.t_inv));
            *exchanged = true;
            *log_det += (1.0 + delta).ln();      
            *av_var = self.make_ti_from_tb().map_err(|e| anyhow!("Failed to make ti from tb: {}", e))?;
            debug_println!("t_inv after make_ti_from_tb: {}", pretty_print!(&self.t_inv));
            debug_println!("t after make_ti_from_tb: {}", pretty_print!(&self.t));
    
            self.transform();
    
    
            debug_println!("block_data.t_x: {}", pretty_print!(&self.t_x));
            debug_println!("block_data.t_block_means: {}", pretty_print!(&self.t_block_means));
            debug_println!("block_data.block_means: {}", pretty_print!(&self.block_means));
        })
    }
    
    fn transform(&mut self) {
        self.t_x = self.x.clone() * self.t_inv.transpose().clone();
        self.t_block_means = self.block_means.clone() * self.t_inv.transpose().clone();
    }

    
}


fn calc_index(i: usize, nc: usize) -> usize {
    i * (nc+1)
}
fn get_range_b(p_mx: &mut Vec<f64>, p_mn: &mut Vec<f64>, vec: &DVector<f64>, k: usize) {
    for i in 0..k {
        p_mx[i] = p_mx[i].max(vec[i]);
        p_mn[i] = p_mn[i].min(vec[i]);
    }
}

/* PermuteB **********************************************************
|	Randomly pemutes the n integers in a[] using the Fike
|	algorithm.  See Fike, "A permutation generation method"  The Computer
|	Journal, 18-1, Feb 75, 21-22.
*/
fn permute_b(a: &mut DVector<u8>, n: u8, random_type: RandomType) -> Result<()> {

    for i in 1..n {
        let rnd = random_type.random();
        //let rnd = 0.5;
        let j = (((1 + i) as f64) * rnd) as i32;
        let temp = a[j as usize];
        a[j as usize] = a[i as usize];
        a[i as usize] = temp;
    }
    
    Ok(())
}




const TOLROT: f64 = 1.0e-12;




/*
fn rmi_from_cmi(column_major_index: usize , width: usize, height: usize) -> usize {
    let row = column_major_index % height;
    let column = column_major_index / height;
    return row * width + column;
}

fn cmi_from_rmi(row_major_index: usize, width: usize, height: usize) -> usize {
    return rmi_from_cmi(row_major_index, height, width);
} */

#[derive(Debug, Clone, Copy)]
enum RandomType {
    Uniform,
    Fixed(f64)
}

impl RandomType {
    fn random(&self) -> f64 {
        match self {
            RandomType::Uniform => rand::random::<f64>(),
            RandomType::Fixed(value) => *value
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct CoincidenceMatrix {
    pub coincidence: DMatrix<usize>,
}

impl CoincidenceMatrix {

    pub fn from_block_array(block_array: &DMatrix<usize>) -> Self {
        //println!("block_array: {}", pretty_print!(&block_array));   
        //let n = block_array.max() + 1;
        let true_n = block_array.max() + 1;
        let n = block_array.nrows();
        let block_size = block_array.ncols();
        let mut coincidence: DMatrix<usize> = DMatrix::zeros(n, n);
        for block_idx in 0..n {
            let block_elements = block_array.row(block_idx);
            for i in 0..block_size {
                let elem_i = block_elements[i];
                // Diagonal counts total appearances
                coincidence[(elem_i, elem_i)] += 1;
                
                // Upper triangle counts pairwise coincidences
                for j in (i+1)..block_size {
                    
                    let elem_j = block_elements[j];
                    debug_println!("i: {}, j: {}, elem_i: {}, elem_j: {}", i, j, elem_i, elem_j);
                    if elem_i < elem_j {
                        coincidence[(elem_i, elem_j)] += 1;
                    } else {
                        coincidence[(elem_j, elem_i)] += 1;
                    }
                }
            }
        }
        Self { coincidence: coincidence.view((0,0),(true_n,true_n)).into() }
    }
}

#[derive(Default, Debug, Clone)]
pub struct BlockArray {
    pub block_array: DMatrix<usize>,
}

impl BlockArray {
    pub fn from_block_array(block_array: &DMatrix<usize>) -> Self {
        Self { block_array: block_array.clone() }
    }

    pub fn as_sorted(&mut self) -> DMatrix<usize> {
        let mut block_array_out = self.block_array.clone(); 
        //println!("block_array before sort: {}", pretty_print!(&block_array_out));
        // first, sort each row ascending
        for mut row in block_array_out.row_iter_mut() {
            let mut swapped = true;

            while swapped {
                swapped = false;
                for i in 0..row.len()-1 {
                    if row[i] > row[i + 1] {
                        row.swap_columns(i, i + 1);
                        swapped = true;
                    }
                }
            }
        }

        for _ in 0..block_array_out.ncols() {
            let mut swapped = true;
            while swapped {
                swapped = false;
                for j in 0..block_array_out.nrows()-1 {
                    let row_j = block_array_out.row(j).clone_owned();
                    let row_j1 = block_array_out.row(j+1).clone_owned();
                    let diff = row_j.cast::<i32>() - row_j1.cast::<i32>();
                    //println!("diff: {}", pretty_print!(&diff));
                    if let Some(first_non_zero) = diff.iter().position(|&x| x != 0) {
                        if diff[first_non_zero] > 0 {
                            block_array_out.swap_rows(j, j+1);
                            swapped = true;
                        }
                    }
                }
            }
        }
        // then, sort each column ascending
        block_array_out
    }
}

#[derive(Debug, Default)]
pub struct BlockResult {
    pub best_log_det: f64,
    pub best_block_array: BlockArray,
    pub best_d: f64,
    pub best_diagonality: f64,
    pub best_coincidence: CoincidenceMatrix
}

// optimize determinant over all blocks using d-criterion





impl BlockDataBuilder {

    fn configure_remaining(&mut self) -> &mut Self {
        let block_data = self;
        if let (Some(n_b), Some(k), Some(block_sizes), Some(n)) = (block_data.n_b, block_data.k, &block_data.block_sizes, block_data.n) {
            block_data.block_means = Some(DMatrix::zeros(n_b as usize, k as usize));
            block_data.t_block_means = Some(DMatrix::zeros(n_b as usize, k as usize));
            block_data.t = Some(DMatrix::zeros(k as usize, k as usize));
            block_data.t_inv = Some(DMatrix::zeros(k as usize, k as usize));
            let n_xb = block_sizes.iter().sum();
            block_data.n_xb = Some(n_xb);
            block_data.rows = Some(DVector::zeros(std::cmp::max(n, n_xb) as usize));
            block_data.b = Some(DMatrix::zeros(n_b as usize, block_sizes[0] as usize));
            block_data.n_t = Some(n);
        }

        block_data
    }

    fn blocksize(&mut self, value: u8) -> &mut Self {
        if let Some(n_b) = self.n_b {
            self.block_sizes = Some(vec![value; n_b as usize]);
            self.max_n = Some(value);
            self
        } else {
            self
        }
    }

    fn v(&mut self, value: u8) -> &mut Self {
        let mut x = DMatrix::zeros(value as usize, value as usize);
        x.fill_diagonal(1.0);
        // remove the first column by subsetting the remaining columns
        let x_sub = x.columns(1, (value - 1) as usize);
        self.x = Some(x_sub.into());
        self.n = Some(value as u8);
        self.k = Some((value - 1) as u8);
        self.t_x =Some(DMatrix::zeros(value as usize, (value - 1) as usize));

        //println!("x: {}", pretty_print!(&x_sub));
        self
    }
} 


pub fn opt_block(v: u8, n_b: u8, block_size: u8,  n_repeats: u8, prohibited_pairs: Vec<(u8, u8)>) -> Result<BlockResult> {
    
    let mut block_data = BlockDataBuilder::default()
        .v(v)
        .n_b(n_b)
        .blocksize(block_size)
        .configure_remaining()
        .prohibited_pairs(prohibited_pairs)
        .build()
        .map_err(|e| anyhow!("Failed to build block_data: {}", e))?;

    
    let block_result = block_data.block_optimize(n_repeats).map_err(|e| anyhow!("Failed to optimize block: {}", e))?;
    Ok(block_result)
}



mod tests {
    use super::*;

    #[allow(unused)]
    fn configure_block_data() -> BlockData {
        let mut block_data = BlockDataBuilder::default()
            .v(7)
            .n_b(7)
            .blocksize(3)
            .random_type(RandomType::Fixed(0.5))
            .configure_remaining()
            .build();
        block_data.unwrap()
    }

    #[allow(unused)]
    fn dm7choose3() -> DMatrix<f64> {
        nalgebra::dmatrix![
            0.0,0.0,0.0,0.0,0.0,0.0;
            1.0,0.0,0.0,0.0,0.0,0.0;
            0.0,1.0,0.0,0.0,0.0,0.0;
            0.0,0.0,1.0,0.0,0.0,0.0;
            0.0,0.0,0.0,1.0,0.0,0.0;
            0.0,0.0,0.0,0.0,1.0,0.0;
            0.0,0.0,0.0,0.0,0.0,1.0;
        ]
    }

    #[test]
    fn test_initialize_b() {
        let mut block_data = configure_block_data();

        block_data.initialize_b().unwrap();

        let expected = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];
        debug_println!("block_data.b: {}", pretty_print!(&block_data.b));
        assert_eq!(block_data.b, expected);
    }

    #[test]
    fn test_form_block_means() {
        //let x = dm7choose3();
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];
        block_data.form_block_means();
        debug_println!("block_data.block_means: {}", pretty_print!(&block_data.block_means));


        let expected = nalgebra::dmatrix![
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0;
            0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0;
            1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0
        ];
        assert_eq!(block_data.block_means, expected);

    }

    #[test]
    fn test_rotate_b() {
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];

        block_data.block_means = nalgebra::dmatrix![    
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0;
            0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0;
            1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0
        ];

        // input 
        let input = nalgebra::dmatrix![
            2.000000, -0.166667, -0.333333,  0.000000,  0.000000, -0.333333;
            0.0,       1.388889, -0.080000, -0.160000, -0.160000, -0.080000;
            0.0,            0.0,  1.768889, -0.010050, -0.198492, -0.695980;
            0.0,            0.0,       0.0,  1.408710, -0.422117, -0.021403;
            0.0,            0.0,       0.0,       0.0,  1.421522, -0.427854;
            0.0,            0.0,       0.0,       0.0,       0.0,  0.651192;
        ];

        block_data.t = input.transpose();
        
        let mut expected = nalgebra::dmatrix![
            2.000000, -0.166667, -0.333333, 0.000000,   0.000000, -0.333333;
            0.0,      1.500000,  -0.074074, -0.296296, -0.074074, -0.074074;
            0.0,           0.0,   1.769547, -0.018605, -0.193023, -0.695349;
            0.0,           0.0,   0.0,       1.756589, -0.465137, -0.031774;
            0.0,           0.0,   0.0,       0.0,       1.434687, -0.421716;
            0.0,           0.0,       0.0,       0.0,       0.0,       0.657029;
        ]; 
        expected.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });

        let vec = nalgebra::dvector![0.000000, -0.333333, 0.000000, 0.666667, -0.333333, 0.000000];

        block_data.rotate_b(&vec, 1.0);
        let out = block_data.t.transpose().apply_into(|x: &mut f64| { *x = (*x * 1000.0).round() });
        assert_eq!(out, expected);
    }


    #[test]
    fn test_reduce_x_to_t() {
        //let x = dm7choose3();
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];

        block_data.block_means = nalgebra::dmatrix![    
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0;
            0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0;
            1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0
        ];
        let (log_det, singular) = block_data.reduce_x_to_t();

        debug_println!("T: {}", pretty_print!(&block_data.t));

        let mut expected_t = nalgebra::dmatrix![
            2.000000, -0.166667, -0.333333, 0.000000,   0.000000, -0.333333;
            0.0,       1.944444, -0.057143, -0.342857, -0.171429, -0.057143;
            0.0,            0.0,  1.771429, -0.021505, -0.198925, -0.693548;
            0.0,            0.0,       0.0,  1.770609, -0.445344, -0.036437;
            0.0,            0.0,       0.0,       0.0,  1.521592, -0.411086;
            0.0,            0.0,       0.0,       0.0,       0.0, 0.659867
        ];
        expected_t.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });
        block_data.t.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });
        assert_eq!((log_det * 1000.0).round(), 2505.0);
        assert_eq!(singular, false);
        assert_eq!(block_data.t, expected_t.transpose());
    }

    
    #[test]
    fn test_make_ti_from_tb() {
        let mut block_data = configure_block_data();
        block_data.t = nalgebra::dmatrix![
            2.000000, -0.166667, -0.333333, 0.000000,   0.000000, -0.333333;
            0.0,       1.944444, -0.057143, -0.342857, -0.171429, -0.057143;
            0.0,            0.0,  1.771429, -0.021505, -0.198925, -0.693548;
            0.0,            0.0,       0.0,  1.770609, -0.445344, -0.036437;
            0.0,            0.0,       0.0,       0.0,  1.521592, -0.411086;
            0.0,            0.0,       0.0,       0.0,       0.0, 0.659867
        ].transpose();

        //block_data.t_inv = block_data.t.transpose();

        let a_var = block_data.make_ti_from_tb().unwrap();

        // Ti is in the lower triangle now

        let mut expected = nalgebra::dmatrix![
            0.707107, 0.0, 0.0, 0.0, 0.0, 0.0;
            0.119523, 0.717137, 0.0, 0.0, 0.0, 0.0;
            0.257603, 0.042934, 0.751343, 0.0, 0.0, 0.0;
            0.048485, 0.258586, 0.016162, 0.751517, 0.0, 0.0;
            0.101746, 0.272416, 0.169029, 0.361033, 0.810683, 0.0;
            0.781205, 0.304621, 0.960265, 0.270228, 0.506063, 1.231039;
        ];
        expected.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });
        block_data.t_inv.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });
        assert_eq!(block_data.t_inv, expected);
        assert_eq!(a_var, 1.0644289011525316);
    }

    #[test]
    fn test_transform() {
        let mut block_data = configure_block_data();
        block_data.t_inv = nalgebra::dmatrix![
            0.707107, 0.0, 0.0, 0.0, 0.0, 0.0;
            0.119523, 0.717137, 0.0, 0.0, 0.0, 0.0;
            0.257603, 0.042934, 0.751343, 0.0, 0.0, 0.0;
            0.048485, 0.258586, 0.016162, 0.751517, 0.0, 0.0;
            0.101746, 0.272416, 0.169029, 0.361033, 0.810683, 0.0;
            0.781205, 0.304621, 0.960265, 0.270228, 0.506063, 1.231039;
        ];
        block_data.block_means = nalgebra::dmatrix![    
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0;
            0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0;
            1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0;
            1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0;
            0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0
        ];
        block_data.t_x = block_data.x.clone() * block_data.t_inv.transpose().clone();
        block_data.t_block_means = block_data.block_means.clone() * block_data.t_inv.transpose().clone();

        let expected_t_x = nalgebra::dmatrix![  
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000;
            0.707107, 0.119523, 0.257603, 0.048485, 0.101746, 0.781205; 
            0.000000, 0.717137, 0.042934, 0.258586, 0.272416, 0.304621; 
            0.000000, 0.000000, 0.751343, 0.016162, 0.169029, 0.960265; 
            0.000000, 0.000000, 0.000000, 0.751517, 0.361033, 0.270228; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.810683, 0.506063; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.231039; 
        ];

        let mut expected_t_block_means = nalgebra::dmatrix![
            0.000000, 0.239046, 0.014311, 0.336701, 0.211149, 0.191616;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.000000, 0.000000, 0.250506, 0.390572, 0.258764; 
            0.000000, 0.000000, 0.250448, 0.005387, 0.326571, 0.899122;
            0.235702, 0.278887, 0.100179, 0.102357, 0.124720, 0.361942;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.239046, 0.014311, 0.336701, 0.481377, 0.360304;
        ];

        expected_t_block_means.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });
        block_data.t_block_means.apply(|x: &mut f64| { *x = (*x * 1000.0).round() });

        debug_println!("block_data.t_x: {}", pretty_print!(&block_data.t_x));
        assert_eq!(block_data.t_x, expected_t_x);
        assert_eq!(block_data.t_block_means, expected_t_block_means);

    }

    #[test]
    fn test_find_delta_block() {
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];
        block_data.t_x = nalgebra::dmatrix![  
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000;
            0.707107, 0.119523, 0.257603, 0.048485, 0.101746, 0.781205; 
            0.000000, 0.717137, 0.042934, 0.258586, 0.272416, 0.304621; 
            0.000000, 0.000000, 0.751343, 0.016162, 0.169029, 0.960265; 
            0.000000, 0.000000, 0.000000, 0.751517, 0.361033, 0.270228; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.810683, 0.506063; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.231039; 
        ];

        block_data.t_inv = nalgebra::dmatrix![
            0.707107, 0.0, 0.0, 0.0, 0.0, 0.0;
            0.119523, 0.717137, 0.0, 0.0, 0.0, 0.0;
            0.257603, 0.042934, 0.751343, 0.0, 0.0, 0.0;
            0.048485, 0.258586, 0.016162, 0.751517, 0.0, 0.0;
            0.101746, 0.272416, 0.169029, 0.361033, 0.810683, 0.0;
            0.781205, 0.304621, 0.960265, 0.270228, 0.506063, 1.231039;
        ];

        block_data.t_block_means = nalgebra::dmatrix![
            0.000000, 0.239046, 0.014311, 0.336701, 0.211149, 0.191616;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.000000, 0.000000, 0.250506, 0.390572, 0.258764; 
            0.000000, 0.000000, 0.250448, 0.005387, 0.326571, 0.899122;
            0.235702, 0.278887, 0.100179, 0.102357, 0.124720, 0.361942;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.239046, 0.014311, 0.336701, 0.481377, 0.360304;
        ];
        let mut new_block = 0;
        let mut x_new = 0;
        debug_println!("block_data.t_x: {}", pretty_print!(&block_data.t_x));
        debug_println!("block_data.t_block_means: {}", pretty_print!(&block_data.t_block_means));
        //find_delta_block(block_data: &mut BlockData, xcur: u8, xnew: &mut u8, cur_block: u8, new_block: &mut u8)  
        block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
        assert_eq!(new_block, 1);
        assert_eq!(x_new, 0);
    }

    #[test]
    fn test_find_delta_block_with_prohibited_pairs() {
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];
        block_data.t_x = nalgebra::dmatrix![  
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000;
            0.707107, 0.119523, 0.257603, 0.048485, 0.101746, 0.781205; 
            0.000000, 0.717137, 0.042934, 0.258586, 0.272416, 0.304621; 
            0.000000, 0.000000, 0.751343, 0.016162, 0.169029, 0.960265; 
            0.000000, 0.000000, 0.000000, 0.751517, 0.361033, 0.270228; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.810683, 0.506063; 
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.231039; 
        ];

        block_data.t_inv = nalgebra::dmatrix![
            0.707107, 0.0, 0.0, 0.0, 0.0, 0.0;
            0.119523, 0.717137, 0.0, 0.0, 0.0, 0.0;
            0.257603, 0.042934, 0.751343, 0.0, 0.0, 0.0;
            0.048485, 0.258586, 0.016162, 0.751517, 0.0, 0.0;
            0.101746, 0.272416, 0.169029, 0.361033, 0.810683, 0.0;
            0.781205, 0.304621, 0.960265, 0.270228, 0.506063, 1.231039;
        ];

        block_data.t_block_means = nalgebra::dmatrix![
            0.000000, 0.239046, 0.014311, 0.336701, 0.211149, 0.191616;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.000000, 0.000000, 0.250506, 0.390572, 0.258764; 
            0.000000, 0.000000, 0.250448, 0.005387, 0.326571, 0.899122;
            0.235702, 0.278887, 0.100179, 0.102357, 0.124720, 0.361942;
            0.235702, 0.039841, 0.336315, 0.021549, 0.090258, 0.990836;
            0.000000, 0.239046, 0.014311, 0.336701, 0.481377, 0.360304;
        ];

        let mut new_block = 0;
        let mut x_new = 0;
        
        // Test Case 1: Prohibited pairs in target block
        block_data.prohibited_pairs = vec![(0, 4)];  // Point 0 can't be with point 4
        let delta1 = block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
        
        if delta1 > DESIGN_TOL {
            // Verify target block doesn't contain prohibited pair
            let target_points: Vec<u8> = block_data.b.row(new_block as usize)
                .iter()
                .take(block_data.block_sizes[new_block as usize] as usize)
                .map(|&x| x as u8)
                .collect();
            
            assert!(!target_points.contains(&4), 
                "Found exchange that would create prohibited pair (0,4) in target block");
        }

        // Test Case 2: Prohibited pairs in current block
        block_data.prohibited_pairs = vec![(2, 5)];  // Point 2 can't be with point 5
        let cur_block = 0;  // Block containing point 2
        let xcur = 1;      // Position of point 2 in block
        
        let delta2 = block_data.find_delta_block(xcur, &mut x_new, cur_block, &mut new_block).unwrap();
        
        if delta2 > DESIGN_TOL {
            // Get the point that would be exchanged into current block
            let incoming_point = block_data.b[(new_block as usize, x_new as usize)] as u8;
            
            // Get remaining points in current block
            let current_block_points: Vec<u8> = block_data.b.row(cur_block as usize)
                .iter()
                .take(block_data.block_sizes[cur_block as usize] as usize)
                .map(|&x| x as u8)
                .filter(|&x| x != 2)  // Exclude point being moved out
                .collect();
            
            // Verify the exchange wouldn't create prohibited pair in current block
            assert!(!current_block_points.contains(&5) || incoming_point != 2,
                "Found exchange that would create prohibited pair (2,5) in current block");
        }

        // Test Case 3: Multiple prohibited pairs
        block_data.prohibited_pairs = vec![(0, 4), (2, 5), (1, 3)];
        let delta3 = block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
        
        if delta3 > DESIGN_TOL {
            // Verify neither block would contain prohibited pairs after exchange
            let target_points: Vec<u8> = block_data.b.row(new_block as usize)
                .iter()
                .take(block_data.block_sizes[new_block as usize] as usize)
                .map(|&x| x as u8)
                .filter(|&x| x != block_data.b[(new_block as usize, x_new as usize)] as u8)
                .collect();
            
            let current_points: Vec<u8> = block_data.b.row(0 as usize)
                .iter()
                .take(block_data.block_sizes[0] as usize)
                .map(|&x| x as u8)
                .filter(|&x| x != 0)  // Exclude point being moved out
                .collect();

            let incoming_point = block_data.b[(new_block as usize, x_new as usize)] as u8;
            
            for &(a, b) in &block_data.prohibited_pairs {
                // Check target block
                assert!(!target_points.contains(&b) || 0 != a, 
                    "Exchange would create prohibited pair in target block");
                
                // Check current block
                assert!(!current_points.contains(&b) || incoming_point != a,
                    "Exchange would create prohibited pair in current block");
            }
        }
    }

    #[test]
    fn test_exchange_blocks() {
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];
        block_data.t = nalgebra::dmatrix![
            2.000000, -0.166667, -0.333333, 0.000000,   0.000000, -0.333333;
            0.0,       1.944444, -0.057143, -0.342857, -0.171429, -0.057143;
            0.0,            0.0,  1.771429, -0.021505, -0.198925, -0.693548;
            0.0,            0.0,       0.0,  1.770609, -0.445344, -0.036437;
            0.0,            0.0,       0.0,       0.0,  1.521592, -0.411086;
            0.0,            0.0,       0.0,       0.0,       0.0, 0.659867
        ];

        let expected_b = nalgebra::dmatrix![
            6,2,4;
            0,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];

        let xcur = 0;
        let xnew = 0;

        let cur_block = 0;
        let mut new_block = 1;


        block_data.exchange_block(xcur, xnew, cur_block, &mut new_block).unwrap();
        debug_println!("block_data.b: {}", pretty_print!(&block_data.b));

        assert_eq!(block_data.b, expected_b);
    }
    
    #[test]
    fn test_block_optimize() {
        let mut block_data = configure_block_data();
        let block_result = block_data.block_optimize(5).unwrap();
        //debug_println!("block_result: {}", pretty_print!(&block_result.best_block_array.cast::<u8>()));

        assert_eq!(block_result.best_log_det, 3.1378770132679095);
    }

    #[test]
    fn test_permute_b() {
        // Create a vector with sequential numbers 0..5
        let mut a = DVector::from_vec(vec![0, 1, 2, 3, 4]);
        let n = 5;
        
        // Use fixed random value of 0.5 to make test deterministic
        let random_type = RandomType::Fixed(0.5);
        
        permute_b(&mut a, n, random_type).unwrap();

        // With fixed random value 0.5, we can predict the exact permutation:
        // i=1: j=floor((1+1)*0.5)=1 -> swap a[1] with a[1] -> [0,1,2,3,4]
        // i=2: j=floor((1+2)*0.5)=1 -> swap a[2] with a[1] -> [0,2,1,3,4]
        // i=3: j=floor((1+3)*0.5)=2 -> swap a[3] with a[2] -> [0,2,3,1,4]
        // i=4: j=floor((1+4)*0.5)=2 -> swap a[4] with a[2] -> [0,2,4,1,3]
        let expected = DVector::from_vec(vec![0, 2, 4, 1, 3]);
        
        assert_eq!(a, expected);
    }

    #[test]
    fn test_no_dup_permute_b() {
        let mut block_data = configure_block_data();
        
        // Initialize rows with sequential numbers
        block_data.rows = DVector::from_vec(vec![6, 5, 4, 3, 2, 1, 0]);
        
        // Set up b matrix with some initial values
        block_data.b = nalgebra::dmatrix![
            0,2,4;
            6,3,1;
            5,0,4;
            3,5,6;
            2,1,0;
            3,6,1;
            5,4,2;
        ];


        // Call no_dup_permute_b with test parameters
        // offset: 0, little_n: 2, bs: 3
        block_data.no_dup_permute_b(0, 2, 3).unwrap();

        // After permutation:
        // 1. The permuted rows should all be different numbers
        // 2. None of the first (bs - little_n) elements should match the values
        //    in block_data.b at the specified offset
        
        // Get the first block row from b matrix
        let first_block = block_data.b.row(0);
        
        // Check that none of the first (bs - little_n = 1) elements in rows
        // match the values in the first block of b
        
        let non_matching = (0..1).all(|j| {
            !first_block.iter().take(2).any(|&val| {
                block_data.rows[j] as i32 == val
            })
        });
 
        assert!(non_matching, "Found matching values between permuted rows and block values");
        
        // Check that all values in rows are unique
        let mut seen = vec![false; 7];
        let all_unique = block_data.rows.iter().all(|&x| {
            if seen[x as usize] {
                false
            } else {
                seen[x as usize] = true;
                true
            }
        }); 

        assert!(all_unique, "Permuted rows contain duplicate values");
    }
}


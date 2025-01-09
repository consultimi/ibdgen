/**
 * opt_block.rs - Design of optimal IBDs with prohibitions
 * Based on the original opt_block.c C code in the AlgDesign R package
 * The original C code is Copyright (c) 2002-2004, Bob Wheeler
 * AlgDesign now maintained by Jerome Braun https://github.com/jvbraun/AlgDesign
 */

use nalgebra::{DMatrix, DVector, Matrix3};
use pretty_print_nalgebra::*;
use anyhow::*;
use derive_builder::Builder;
use crate::block_array::BlockArray;

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

    /* n is the number of rows in x */
    n: u8,

    /* k is the number of columns in x */
    k: u8, 

    /* n_t is the number of rows to use. If init_rows is true, n_t = n_xb, otherwise n_t = n */
    n_t: usize,

    /* n_xb is the number of rows in x that are used in the blocks */
    n_xb: usize,           

    /* n_b is the number of blocks */
    n_b: usize,           
    
    #[builder(default = "RandomType::Uniform")]
    random_type: RandomType, // random_type is the type of randomization
    
    block_size: usize, // block_size is the size of each block

    /* rows is a vector of row indices */
    rows: DVector<usize>, 

    /* prohibited_pairs is a vector of prohibited pairs */
    #[builder(default = "vec![]")]
    prohibited_pairs: Vec<(usize, usize)>,

    /* moments contain the calculation for potential improvements */
    moments: Matrix3<f64>,

    /* A matrix used to store the diffs used in find_delta_block */
    diffs: DMatrix<f64>,
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

    fn initialize_block_array(&mut self, block_array: &mut Vec<usize>)  {

        for i in 0..self.n {
            self.rows[i as usize] = i as usize;
        }
    
        let mut l = 0;
        let mut m = 0;
        let n_t = self.n;
        for _ in 0..self.n_b {
            for _ in 0..self.block_size {
                if l >= n_t {
                    l = 0;
                }
                block_array[m as usize] = self.rows[l as usize] + 1;
                m += 1;
                l += 1;
            }
        }
    }

    /*
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
    } */
    
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
            block.copy_from(&(DMatrix::from_rows(&out).row_sum() / self.block_size as f64));
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
        for i in 0..self.n_b * self.block_size {
            self.b[((i / self.block_size) as usize, (i % self.block_size) as usize)] = -1;
        }

        let mut l = 0;  // Index into rows array
        
        // 4. For each block
        for i in 0..self.n_b {
            let mut retry_count = 0;
            
            // 5. Fill each position in the block
            let mut j = 0;
            while j < self.block_size {
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
                    let existing_point = self.b[(i as usize, k as usize)] as usize;
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

    
    fn find_delta_block(&mut self, xcur: usize, xnew: &mut usize, cur_block: usize, new_block: &mut usize) -> Result<f64> {
        const DELTA_TOL: f64 = 1e-12;  // Minimum improvement threshold
        let mut delta = 0.0;  // Tracks best improvement found
        //let mut g_vec: SVector<f64, 3> = SVector::from_vec(vec![0.0, 1.0, 0.0]);
        //let mut mi_vec: SVector<f64, 3> = SVector::from_vec(vec![0.0, 0.0, 0.0]);
        let cur_treatment_rowno = self.b[(cur_block as usize, xcur as usize)] as usize;

        let fi = self.t_x.row(cur_treatment_rowno);
        let fmi = self.t_block_means.row(cur_block as usize);

        // Loop through all blocks except current
        for i in 0..self.n_b {
            if i == cur_block {
                continue;
            }
            let fmj = self.t_block_means.row(i as usize);
            let fmj_fmi_diff = fmj - fmi;
            //let vectors = DMatrix::from_rows(&[fmj_fmi_diff, fmj_fmi_diff]).transpose();
            self.diffs.set_column(0, &fmj_fmi_diff.transpose());
                // Try exchanging with each point in candidate block
            for j in 0..self.block_size {
                let candidate_treatment_rowno = self.b[(i as usize, j as usize)] as usize;

                // Skip if prohibited pairs exist
                if !self.prohibited_pairs.is_empty() {
                    let is_valid = self.check_prohibited_pairs(xcur, cur_block, cur_treatment_rowno, i, j, candidate_treatment_rowno);

                    if !is_valid {
                        continue;
                    }
                }

                let fj = self.t_x.row(candidate_treatment_rowno as usize);
                let fj_fi_diff = fj - fi;
                self.diffs.set_column(1, &fj_fi_diff.transpose());
                //self.moments[1] = fmj_fmi_diff.component_mul(&fmj_fmi_diff).sum();
                //println!("fmj_fmi_diff: {}", pretty_print!(&fmj_fmi_diff));
                //println!("fj_fi_diff: {}", pretty_print!(&fj_fi_diff));

                
                //println!("vectors: {}", &vectors);
                // Multiply with its transpose to get a 2×2 matrix of dot products
                let dot_products = self.diffs.tr_mul(&self.diffs);
                //let dot_products = vectors.clone() * vectors.transpose();

                //println!("dot_products: {}", dot_products);
                self.moments[1] = dot_products[0]; 
                self.moments[4] = dot_products[2]; 
                self.moments[7] = dot_products[3]; 

                // Combine geometric and moment terms
                
                
                // Combine geometric and moment terms
                self.moments.set_row(2, &(self.moments.row(0) + self.moments.row(1)));
                
                let d = -(1.0 + self.moments[2] * self.moments[8] - self.moments[5] * self.moments[5]);

                if (d - delta) > DELTA_TOL {
                    delta = d;
                    *new_block = i;
                    *xnew = j;
                }
            }
        }

        Ok(delta)
    }

    fn check_prohibited_pairs(&self, xcur: usize, cur_block: usize, cur_treatment_rowno: usize, i: usize, j: usize, candidate_treatment_rowno: usize) -> bool {
        // Check target block (only need to check up to j)
        for k in 0..self.block_size {
            if k != j {  // Skip the point being exchanged
                let point = self.b[(i as usize, k as usize)] as usize;
                if self.prohibited_pairs.iter().any(|&(a, b)| 
                    (cur_treatment_rowno == a && point == b) || 
                    (cur_treatment_rowno == b && point == a)
                ) {
                    return false;
                }
            }
        }
    
        for k in 0..self.block_size {
            if k != xcur {  // Skip the point being exchanged
                let point = self.b[(cur_block as usize, k as usize)] as usize;
                if self.prohibited_pairs.iter().any(|&(a, b)| 
                    (candidate_treatment_rowno == a && point == b) || 
                    (candidate_treatment_rowno == b && point == a)
                ) {
                    return false;
                }
            }
        }
        true
    }

    
    fn exchange_block(&mut self, xcur: usize, xnew: usize, cur_block: usize, new_block: &mut usize) -> Result<()> {
        let row_no_i = self.b[(cur_block as usize, xcur as usize)] as usize;
        //let ni = self.block_size;

        let x_clone = self.x.clone();
        let xri = x_clone.row(row_no_i);
        let xmi = self.block_means.row(cur_block as usize);
        //debug_println!("xri: {}\nxmi: {}\nrowNoi: {}", pretty_print!(&xri), pretty_print!(&xmi), row_no_i);

        let row_no_j = self.b[(*new_block as usize, xnew as usize)] as usize;
        // Handle normal block exchange case

        let xrj = x_clone.row(row_no_j);
        let xmj = self.block_means.row(*new_block as usize);
        debug_println!("xmi: {}\nxmj: {}\nxri: {}\nxrj: {}\nrowNoj: {}", pretty_print!(&xmi), pretty_print!(&xmj), pretty_print!(&xri), pretty_print!(&xrj), row_no_j);
        //let nj = self.block_size;
        //let c = (ni + nj) as f64 / (ni * nj) as f64;

        // vec = xmj - xmi
        let mut vec = (xmj - xmi).transpose();
        self.rotate_b(&vec, 1.0);
        vec -= (xrj - xri).transpose();
        self.rotate_b(&vec, -1.0);
        vec = (xrj - xri).transpose();
        self.rotate_b(&vec, 1.0 - self.moments[0]);
        
        // Update block means
        for i in 0..self.k {
            //let idx = cmi_from_rmi((cur_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            //self.block_means[idx as usize] += newsum;
            self.block_means[(cur_block as usize, i as usize)] += (xrj[i as usize] - xri[i as usize]) / self.block_size as f64;
            //let idx = cmi_from_rmi((*new_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            self.block_means[(*new_block as usize, i as usize)] += (xri[i as usize] - xrj[i as usize]) / self.block_size as f64;
        }

        //println!("new_block: {}, xnew: {}, b before exchange: {}", *new_block, xnew, pretty_print!(&self.b));

        self.b[(*new_block as usize, xnew as usize)] = row_no_i as i32;
        self.b[(cur_block as usize, xcur as usize)] = row_no_j as i32;

        //println!("cur_block: {}, xcur: {}, b after exchange: {}", cur_block, xcur, pretty_print!(&self.b));

        Ok(())
    }

    fn block_optimize(&mut self, n_repeats: usize) -> Result<BlockResult> {

        let mut block_array: Vec<usize> = vec![0; self.n_b as usize * self.block_size as usize];
        let mut best_log_det = 0.0;
        let mut best_block_array = DMatrix::zeros(self.n_b as usize, self.block_size as usize);
        let mut xnew = 0;
        let mut new_block = 0;
        let mut av_var = 0.0;
        // b is a matrix of block factors. ncols is max(blocksizes)
        self.initialize_block_array(&mut block_array);
    
        for repeat_num in 0..n_repeats {
            debug_println!("REPEAT NUMBER: {}", repeat_num + 1);

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
                        for xcur in 0..self.block_size {
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

    fn try_exchange_block(&mut self, xnew: &mut usize, new_block: &mut usize, av_var: &mut f64, log_det: &mut f64, exchanged: &mut bool, cur_block: usize, xcur: usize) -> Result<(), Error> {
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
fn permute_b(a: &mut DVector<usize>, n: usize, random_type: RandomType) -> Result<()> {

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


    pub fn r(&self) -> f64 {
        self.coincidence.diagonal().cast::<f64>().mean()
    }

    pub fn variance(&self) -> f64 {
        let mut upper_try_f: nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self.coincidence.upper_triangle().cast::<f64>();
        upper_try_f.fill_lower_triangle_with_upper_triangle();
        upper_try_f.fill_diagonal(self.lambda());
        upper_try_f.variance()
    }

    pub fn lambda(&self) -> f64 {
        let mut upper_try_f: nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self.coincidence.upper_triangle().cast::<f64>();
        upper_try_f.fill_diagonal(0.0);

        let cells = upper_try_f.ncols() * (upper_try_f.ncols() - 1) / 2;
        let lambda = upper_try_f.sum() as f64 / (cells as f64);
        lambda
    }

    pub fn is_bibd(&self) -> bool {
        let diagonal_f = self.coincidence.diagonal().cast::<f64>();
        (diagonal_f.max() == diagonal_f.min()) && (self.variance() == 0.0)
    }

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
        if let (Some(n_b), Some(k), Some(block_size), Some(n)) = (block_data.n_b, block_data.k, &block_data.block_size, block_data.n) {
            block_data.block_means = Some(DMatrix::zeros(n_b as usize, k as usize));
            block_data.t_block_means = Some(DMatrix::zeros(n_b as usize, k as usize));
            block_data.t = Some(DMatrix::zeros(k as usize, k as usize));
            block_data.t_inv = Some(DMatrix::zeros(k as usize, k as usize));
            let n_xb = block_size * n_b;
            block_data.n_xb = Some(n_xb);
            block_data.rows = Some(DVector::zeros(std::cmp::max(n as usize, n_xb as usize)));
            block_data.b = Some(DMatrix::zeros(n_b as usize, *block_size as usize));
            block_data.n_t = Some(n as usize);
            block_data.moments = Some(Matrix3::new(
                (2 * block_size) as f64 / (block_size * block_size) as f64, 1.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            ));
            block_data.diffs = Some(DMatrix::zeros(k as usize, 2));
        }

        block_data
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


pub fn opt_block(v: u8, n_b: usize, block_size: usize,  n_repeats: usize, prohibited_pairs: Vec<(usize, usize)>) -> Result<BlockResult> {
    
    let mut block_data = BlockDataBuilder::default()
        .v(v)
        .n_b(n_b)
        .block_size(block_size)
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
            .block_size(3)
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

        // Test Case 1: Single prohibited pair in target block
        {
            block_data.prohibited_pairs = vec![(0, 4)];
            let delta = block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
            
            if delta > DESIGN_TOL {
                // Verify the proposed exchange doesn't put 0 and 4 together
                let target_block = block_data.b.row(new_block as usize);
                let has_four = (0..block_data.block_size)
                    .any(|k| k != x_new && target_block[k as usize] == 4);
                assert!(!has_four, "Exchange would create prohibited pair (0,4) in target block");
            }
        }

        // Test Case 2: Single prohibited pair in current block
        {
            block_data.prohibited_pairs = vec![(2, 5)];
            let cur_block = 0;
            let xcur = 1;  // Position of point 2
            
            let delta = block_data.find_delta_block(xcur, &mut x_new, cur_block, &mut new_block).unwrap();
            
            if delta > DESIGN_TOL {
                let incoming_point = block_data.b[(new_block as usize, x_new as usize)] as u8;
                let has_five = (0..block_data.block_size)
                    .any(|k| k != xcur && block_data.b[(cur_block as usize, k as usize)] == 5);
                
                assert!(!has_five || incoming_point != 2, 
                    "Exchange would create prohibited pair (2,5) in current block");
            }
        }

        // Test Case 3: Multiple prohibited pairs affecting both blocks
        {
            block_data.prohibited_pairs = vec![(0, 4), (2, 5), (1, 3)];
            let delta = block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
            
            if delta > DESIGN_TOL {
                let incoming_point = block_data.b[(new_block as usize, x_new as usize)] as usize;
                
                // Check target block
                for k in 0..block_data.block_size {
                    if k != x_new {
                        let point = block_data.b[(new_block as usize, k as usize)] as usize;
                        assert!(!block_data.prohibited_pairs.iter().any(|&(a, b)| 
                            (0 == a && point == b) || (0 == b && point == a)
                        ), "Exchange would create prohibited pair in target block");
                    }
                }
                
                // Check current block
                for k in 0..block_data.block_size {
                    if k != 0 {
                        let point = block_data.b[(0, k as usize)] as usize;
                        assert!(!block_data.prohibited_pairs.iter().any(|&(a, b)| 
                            (incoming_point == a && point == b) || 
                            (incoming_point == b && point == a)
                        ), "Exchange would create prohibited pair in current block");
                    }
                }
            }
        }

        // Test Case 4: Performance - many prohibited pairs
        {
            // Create many prohibited pairs
            block_data.prohibited_pairs = (0..6)
                .flat_map(|i| (i+1..7).map(move |j| (i, j)))
                .collect();
            
            use std::time::Instant;
            let start = Instant::now();
            let delta = block_data.find_delta_block(0, &mut x_new, 0, &mut new_block).unwrap();
            let duration = start.elapsed();
            
            // Should complete in reasonable time (e.g., < 1ms)
            assert!(duration.as_micros() < 1000, 
                "find_delta_block took too long with many prohibited pairs");
            
            // Verify result if exchange found
            if delta > DESIGN_TOL {
                let incoming_point = block_data.b[(new_block as usize, x_new as usize)] as usize;
                
                // Quick check of one prohibited pair in each block
                let target_has_violation = (0..block_data.block_size)
                    .any(|k| k != x_new && block_data.prohibited_pairs.contains(&(
                        0,
                        block_data.b[(new_block as usize, k as usize)] as usize
                    )));
                    
                let current_has_violation = (0..block_data.block_size)
                    .any(|k| k != 0 && block_data.prohibited_pairs.contains(&(
                        incoming_point,
                        block_data.b[(0, k as usize)] as usize
                    )));
                    
                assert!(!target_has_violation && !current_has_violation,
                    "Exchange with many prohibited pairs created violation");
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

    /*
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
     */
}


use nalgebra::{DMatrix, DVector};
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
            n_xb: 0,
            n_b: 0,
            block_sizes: vec![0; 0],
            n_repeat_counts: 0,
            rows: DVector::zeros(0),
        }
    }



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
                debug_println!("i: {}, j: {}, k_index: {}, s: {}, t_vec[{}]: {}, c: {}, r: {}", i, j, k_index, s, j, t_vec[j as usize], c, self.t[k_index]);
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
        //debug_println!("no_dup_permute_b called with little_n: {}, bs: {}, block_data.rows: {}, offset: {}", little_n, bs, &block_data.rows, offset);
        //debug_println!("b: {}", pretty_print!(&block_data.b));
        loop {
            //dbg!(&block_data.rows);
            let mut nodup = true;
            permute_b(&mut self.rows, self.n).map_err(|e| anyhow!("Failed to permute rows: {}", e))?;
            for i in 0..little_n {
                debug_println!("bs: {}, little_n: {}, offset: {}, i: {}", bs, little_n, offset, i);
                let index = offset * self.max_n + i;
                let cur_val = self.b[cmi_from_rmi(index as usize, self.max_n as usize, self.n_b as usize) as usize];
                //let index = offset * block_data.k + i;
                //debug_println!("index: {}", index);
                //let cur_val = block_data.b[index as usize];
                for j in 0..(bs - little_n) {
                    if self.rows[j as usize] as f64 == cur_val {
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

        for i in 0..self.n_t {
            self.rows[i as usize] = i;
        }


        permute_b(&mut self.rows, self.n_t).map_err(|e| anyhow!("Failed to permute rows: {}", e))?;

        /*	for (i=0;i<nB*MAXN;i++)
            B[i]=-1; */
        for i in 0..self.n_b * self.max_n {
            self.b[i as usize] = -1.0;
        }

        let mut l = 0;
        for i in 0..self.n_b {
            let bs = self.block_sizes[i as usize];
            debug_println!("bs: {}", bs);
            for j in 0..bs {
                if l >= self.n_t {
                    l = 0;
                    self.no_dup_permute_b(i, j, bs)?;
                }
                let index = (i * self.max_n + j) as usize;
                let column_major_index = cmi_from_rmi(index, self.max_n as usize, self.n_b as usize);
                //debug_println!("column_major_index: {}, block_data.rows[l as usize]: {}", column_major_index, block_data.rows[l as usize]);
                self.b[column_major_index] = self.rows[l as usize] as f64;
                l += 1;
            }
        }
        debug_println!("block_data.b after initialize_b: {}", pretty_print!(&self.b));
        Ok(())
    }

    
    fn find_delta_block(&mut self, xcur: u8, xnew: &mut u8, cur_block: u8, new_block: &mut u8) -> Result<f64> {
        const DELTA_TOL: f64 = 1e-12;  // Minimum improvement threshold
        let mut delta = 0.0;  // Tracks best improvement found
        
        debug_println!("find_delta_block called with xcur: {}, xnew: {}, cur_block: {}, new_block: {}", xcur, xnew, cur_block, new_block);
        // Get current point's row number and block size
        //let cur_row_no = block_data.b[(cur_block * block_data.max_n + xcur) as usize] as usize;
        let cur_row_no = self.b[cmi_from_rmi((cur_block * self.max_n + xcur) as usize, self.max_n as usize, self.n_b as usize) as usize] as usize;
        debug_println!("cur_row_no: {}", cur_row_no);
        let ni = self.block_sizes[cur_block as usize];
        //debug_println!("transposed block means: {}", &block_data.t_block_means);

        // Get pointers to current point and its block mean
        let fi = self.t_x.row(cur_row_no);
        let fmi = self.t_block_means.row(cur_block as usize);
        let b_transpose = self.b.transpose();
        debug_println!("t_block_means inside find_delta_blocks: {}", pretty_print!(&self.t_block_means));
        debug_println!("b_transpose inside find_delta_blocks: {}", pretty_print!(&b_transpose));
        debug_println!("b inside find_delta_blocks: {}", pretty_print!(&self.b));
        // Loop through all blocks except current
        for i in 0..self.n_b {
            if i != cur_block {
                let nj = self.block_sizes[i as usize];
                
                // Calculate geometric coefficient based on block sizes
                let gi0 = (ni + nj) as f64 / (ni * nj) as f64;
                let gi1 = 1.0;
                let gi2 = 0.0;

                let fmj = self.t_block_means.row(i as usize);
                
                // Calculate squared distance between block means
                let mut g = 0.0;
                for l in 0..self.k {
                    debug_println!("fmj[{}]: {}, fmi[{}]: {}", l, fmj[l as usize], l, fmi[l as usize]);
                    let dif = fmj[l as usize] - fmi[l as usize];
                    g += dif * dif;
                }
                let mi0 = g;

                
                // Try exchanging with each point in candidate block
                for j in 0..nj {
                    let row_no = b_transpose[(i * self.max_n + j) as usize] as usize;
                    let fj = self.t_x.row(row_no);
                    
                    let mut g = 0.0;
                    let mut h = 0.0;
                    
                    // Calculate cross terms between means and points
                    for l in 0..self.k {
                        //debug_println!("row_no: {}, i: {}, j: {}, l: {}, fmi: {}, fmj: {}, fi: {}, fj: {}", row_no, i, j, l, fmi[l as usize], fmj[l as usize], fi[l as usize], fj[l as usize]);

                        let dif1 = fmj[l as usize] - fmi[l as usize];
                        let dif2 = fj[l as usize] - fi[l as usize];
                        g += dif1 * dif2;
                        h += dif2 * dif2;
                    }
                    let mi1 = g;
                    let mi2 = h;
                    
                    // Combine geometric and moment terms
                    let m1i0 = gi0 + mi0;
                    let m1i1 = gi1 + mi1;
                    let m1i2 = gi2 + mi2;

                    // Calculate improvement in criterion
                    let d = -(1.0 + m1i0 * m1i2 - m1i1 * m1i1);
                    debug_println!("d: {}, i: {}, j: {}", d, i, j);
                    // Update best exchange if improvement is large enough
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
        let mut vec = DVector::zeros(self.k as usize);
        //let b_transpose = block_data.b.transpose();
        let row_no_i = self.b[cmi_from_rmi(
            (cur_block * self.max_n + xcur) as usize,
            self.max_n as usize,
            self.n_b as usize
        ) as usize] as usize;
        let ni = self.block_sizes[cur_block as usize];

        let x_clone = self.x.clone();
        let xri = x_clone.row(row_no_i);
        let xmi = self.block_means.row(cur_block as usize);
        //debug_println!("xri: {}\nxmi: {}\nrowNoi: {}", pretty_print!(&xri), pretty_print!(&xmi), row_no_i);
        // Handle normal block exchange case
        let row_no_j = self.b[cmi_from_rmi(
            (*new_block * self.max_n + xnew) as usize,
            self.max_n as usize,
            self.n_b as usize
        ) as usize] as usize;
        let xrj = x_clone.row(row_no_j);
        let xmj = self.block_means.row(*new_block as usize);
        debug_println!("xmi: {}\nxmj: {}\nxri: {}\nxrj: {}\nrowNoj: {}", pretty_print!(&xmi), pretty_print!(&xmj), pretty_print!(&xri), pretty_print!(&xrj), row_no_j);
        let nj = self.block_sizes[*new_block as usize];
        let c = (ni + nj) as f64 / (ni * nj) as f64;

        // vec = xmj - xmi
        for i in 0..self.k {
            vec[i as usize] = xmj[i as usize] - xmi[i as usize];
        }
        debug_println!("t before rotate: {}", pretty_print!(&self.t));
        self.rotate_b(&vec, 1.0);
        debug_println!("t after rotate: {}", pretty_print!(&self.t));

        // vec -= xrj - xri
        for i in 0..self.k {
            vec[i as usize] -= xrj[i as usize] - xri[i as usize];
        }

        debug_println!("t before rotate: {}", pretty_print!(&self.t));
        self.rotate_b(&vec, -1.0);
        debug_println!("t after rotate: {}", pretty_print!(&self.t));
        // vec = xrj - xri
        for i in 0..self.k {
            vec[i as usize] = xrj[i as usize] - xri[i as usize];
        }
        debug_println!("t before rotate: {}", pretty_print!(&self.t));
        self.rotate_b(&vec, 1.0 - c);
        debug_println!("t after rotate: {}", pretty_print!(&self.t));

        debug_println!("block_data.block_means before update: {}", pretty_print!(&self.block_means));
        // Update block means
        debug_println!("block_data.k: {}", self.k);
        for i in 0..self.k {
            let idx = cmi_from_rmi((cur_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            let newsum = (xrj[i as usize] - xri[i as usize]) / ni as f64;
            debug_println!("idx: {}, xrj[i]: {}, xri[i]: {}, ni: {}, newsum: {}", idx, xrj[i as usize], xri[i as usize], ni, newsum);
            self.block_means[idx as usize] += newsum;
            let idx = cmi_from_rmi((*new_block as usize * self.k as usize + i as usize) as usize, self.k as usize, self.n_b as usize);
            //let idx = row_no_j as u8 + i as u8;
            let newsum = (xri[i as usize] - xrj[i as usize]) / nj as f64;
            debug_println!("idx: {}, xrj[i]: {}, xri[i]: {}, nj: {}, newsum: {}", idx, xrj[i as usize], xri[i as usize], nj, newsum);
            self.block_means[idx as usize] = self.block_means[idx as usize] + newsum;
        }

        debug_println!("block_data.block_means after update: {}", pretty_print!(&self.block_means));

        self.b[cmi_from_rmi(
            (*new_block * self.max_n + xnew) as usize,
            self.max_n as usize,
            self.n_b as usize
        ) as usize] = row_no_i as f64;

        self.b[cmi_from_rmi(
            (cur_block * self.max_n + xcur) as usize,
            self.max_n as usize,
            self.n_b as usize
        ) as usize] = row_no_j as f64;

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
            debug_println!("t_inv after beginning of loop: {}", pretty_print!(&self.t_inv));
            debug_println!("t after first beginning of loop: {}", pretty_print!(&self.t));
            debug_println!("block_data.t_x: {}", pretty_print!(&self.t_x));
            debug_println!("block_data.t_block_means: {}", pretty_print!(&self.t_block_means));
            debug_println!("block_data.block_means: {}", pretty_print!(&self.block_means));
            debug_println!("block_data.b: {}", pretty_print!(&self.b));
    
            self.initialize_b().map_err(|e| anyhow!("Failed to initialize b: {}", e))?;
    
            debug_println!("block_data.b after initialize_b: {}", pretty_print!(&self.b));
    
            //dbg!(&block_data.b);
            self.form_block_means();
            let (mut log_det, singular) = self.reduce_x_to_t();
            if singular {
                return Err(anyhow!("Singular matrix"));
            } else {
                debug_println!("t_inv before first make_ti_from_tb: {}", pretty_print!(&self.t_inv));
                av_var = self.make_ti_from_tb().map_err(|e| anyhow!("Failed to make ti from tb: {}", e))?;
                debug_println!("t_inv after first make_ti_from_tb: {}", pretty_print!(&self.t_inv));
                debug_println!("t after first make_ti_from_tb: {}", pretty_print!(&self.t));
                self.t_x = self.x.clone() * self.t_inv.transpose().clone();
                self.t_block_means = self.block_means.clone() * self.t_inv.transpose().clone();
    
                //debug_println!("420 block_data.t_x: {}", pretty_print!(&block_data.t_x));
                //debug_println!("421block_data.t_block_means: {}", pretty_print!(&block_data.t_block_means));
                loop {  
                    let mut exchanged = false;
                    let mut cur_block = 0;
                    loop {
                        for xcur in 0..self.block_sizes[cur_block as usize] {
                            debug_println!("BEING LOOP xcur: {}, curBlock: {}, newBlock: {}", xcur, cur_block, new_block);
                            let delta = self.find_delta_block(xcur, &mut xnew, cur_block, &mut new_block).map_err(|e| anyhow!("Failed to find delta block: {}", e))?;
                            debug_println!("delta: {}", delta);
                            if delta < 10.0 && delta > DESIGN_TOL {
    
                                debug_println!("t before exchange: {}", pretty_print!(&self.t));
                                debug_println!("t_inv before exchange: {}", pretty_print!(&self.t_inv));
                                self.exchange_block(xcur, xnew, cur_block, &mut new_block).map_err(|e| anyhow!("Failed to exchange block: {}", e))?;
                                debug_println!("t_inv after exchange: {}", pretty_print!(&self.t_inv));
                                exchanged = true;
                                log_det += (1.0 + delta).ln();      
                                av_var = self.make_ti_from_tb().map_err(|e| anyhow!("Failed to make ti from tb: {}", e))?;
                                debug_println!("t_inv after make_ti_from_tb: {}", pretty_print!(&self.t_inv));
                                debug_println!("t after make_ti_from_tb: {}", pretty_print!(&self.t));
    
                                //transform(block_data)?;
                                self.t_x = self.x.clone() * self.t_inv.transpose().clone();
                                self.t_block_means = self.block_means.clone() * self.t_inv.transpose().clone();
    
                                debug_println!("block_data.t_x: {}", pretty_print!(&self.t_x));
                                debug_println!("block_data.t_block_means: {}", pretty_print!(&self.t_block_means));
                                debug_println!("block_data.block_means: {}", pretty_print!(&self.block_means));
                            }
                        }
                        if cur_block == self.n_b - 1 {
                            break;
                        } else {
                            cur_block += 1;
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
        println!("block_data.k: {}", self.k);
        println!("block_data.n_xb: {}", self.n_xb);
        // 	*D=exp(logDbest/(double)k)/(double)Nxb;
        let best_d = (best_log_det / self.k as f64).exp() / self.n_xb as f64;
        let best_diagonality = 1.0 / (best_d * av_var * self.n_xb as f64);
        debug_println!("best_log_det: {}", best_log_det);
        debug_println!("best_block_array: {}", best_block_array);
        Ok(BlockResult { best_log_det, best_block_array, best_d, best_diagonality })
    }

    fn create_coincidence_matrix(&self, block_result: &BlockResult) -> DMatrix<usize> {
        let mut coincidence: DMatrix<usize> = DMatrix::zeros(self.n as usize, self.n as usize);
        for block_idx in 0..self.n_b as usize {
            let block_size = self.block_sizes[block_idx] as usize;
            let block_elements = block_result.best_block_array.row(block_idx);
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
    
        coincidence
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
fn permute_b(a: &mut DVector<u8>, n: u8) -> Result<()> {

    for i in 1..n {
        let rnd = rand::random::<f64>();
        //let rnd = 0.5;
        let j = (((1 + i) as f64) * rnd) as i32;
        let temp = a[j as usize];
        a[j as usize] = a[i as usize];
        a[i as usize] = temp;
    }
    
    Ok(())
}




const TOLROT: f64 = 1.0e-12;





fn rmi_from_cmi(column_major_index: usize , width: usize, height: usize) -> usize {
    let row = column_major_index % height;
    let column = column_major_index / height;
    return row * width + column;
}

fn cmi_from_rmi(row_major_index: usize, width: usize, height: usize) -> usize {
    return rmi_from_cmi(row_major_index, height, width);
}





#[derive(Debug)]
struct BlockResult {
    best_log_det: f64,
    best_block_array: DMatrix<usize>,
    best_d: f64,
    best_diagonality: f64
}

// optimize determinant over all blocks using d-criterion


#[derive(Builder, Debug)]
#[builder(build_fn(error = "anyhow::Error"))]
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
    n_xb: u8,           // n_xb is the number of rows in x that are used in the blocks
    n_b: u8,            // n_b is the number of blocks
    
    #[builder(setter(custom))]
    block_sizes: Vec<u8>, // block_sizes is a vector of block sizes
    n_repeat_counts: u8, // n_repeat_counts is the number of repeats
    rows: DVector<u8>,  // rows is a vector of row indices
}


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
} 


pub fn opt_block(x_i: DMatrix<f64>, n_b: u8, block_size: u8,  n_repeats: u8) -> Result<()> {
    
    let mut block_data = BlockDataBuilder::default()
        .x(x_i.clone())
        .n(x_i.nrows() as u8)
        .k(x_i.ncols() as u8)
        .n_b(n_b)
        .blocksize(block_size)
        .n_repeat_counts(n_repeats)
        .configure_remaining()
        .build()
        .map_err(|e| anyhow!("Failed to build block_data: {}", e))?;

    
    let block_result = block_data.block_optimize(n_repeats).map_err(|e| anyhow!("Failed to optimize block: {}", e))?;
    println!("block_result: {}", pretty_print!(&block_result.best_block_array.add_scalar(1)));
    println!("block_result.best_log_det: {}", block_result.best_log_det);
    println!("block_result.best_d: {}", block_result.best_d);
    println!("block_result.best_diagonality: {}", block_result.best_diagonality);
    // Create coincidence matrix to store pairwise counts
    
    let coincidence = block_data.create_coincidence_matrix(&block_result);
    //debug_println!("coincidence: {}", pretty_print!(&coincidence));
    // For each block, count coincidences between all pairs of elements
    

    println!("Coincidence matrix:");
    println!("{}", coincidence);
    Ok(())
}



mod tests {
    use super::*;

    #[allow(unused)]
    fn configure_block_data() -> BlockData {
        let mut block_data = BlockData::new();
        let x_i = dm7choose3();
        let n_b = 7;
        let block_sizes = vec![3, 3, 3, 3, 3, 3, 3];
        block_data.x = x_i.clone();
        block_data.n = x_i.nrows() as u8;
        block_data.k = x_i.ncols() as u8;
        block_data.block_sizes = block_sizes.clone();
        block_data.block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
        block_data.t_block_means = DMatrix::zeros(n_b as usize, block_data.k as usize);
        block_data.t = DMatrix::zeros(block_data.k as usize, block_data.k as usize);
        block_data.t_inv = DMatrix::zeros(block_data.k as usize, block_data.k as usize);
        block_data.n_repeat_counts = 1;
        block_data.max_n = *block_sizes.iter().max().unwrap();
        block_data.n_xb = block_data.block_sizes.iter().sum();
        block_data.n_b = n_b;
        block_data.rows = DVector::zeros(block_data.n as usize);
        block_data.b = DMatrix::zeros(n_b as usize, block_data.max_n as usize);
        block_data.n_t = block_data.n;
        block_data
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
        let x = dm7choose3();
        let mut block_data = configure_block_data();
        let mut block_array = vec![0; block_data.n_b as usize * block_data.max_n as usize];

       
        block_data.initialize_block_array(&mut block_array);
        //assert_eq!(block_array, expected.to_vec());

        block_data.initialize_b().unwrap();

        let expected = nalgebra::dmatrix![
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
        ];
        debug_println!("block_data.b: {}", pretty_print!(&block_data.b));
        assert_eq!(block_data.b, expected);
    }

    #[test]
    fn test_form_block_means() {
        //let x = dm7choose3();
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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
    fn test_exchange_blocks() {
        let mut block_data = configure_block_data();
        block_data.b = nalgebra::dmatrix![
            0.0,2.0,4.0;
            6.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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
            6.0,2.0,4.0;
            0.0,3.0,1.0;
            5.0,0.0,4.0;
            3.0,5.0,6.0;
            2.0,1.0,0.0;
            3.0,6.0,1.0;
            5.0,4.0,2.0;
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

}


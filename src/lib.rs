mod ibd;
mod block_array;
mod random_type;
mod coincidence_matrix;
mod block_result;

use ibd::*;
use block_result::BlockResult;
use anyhow::Error;

pub fn find_best_ibd(v: u8, n_b: usize, block_size: usize, n_repeats: usize, iter: usize, prohibited_pairs: Vec<(usize, usize)>) -> Result<(BlockResult, usize, f64), Error> {
    let mut min_d = 0.0;
    let mut best_solution = BlockResult::default();
    let mut best_iter = 0;
    for i in 1..=iter {
        let result = ibdgen(
            v, 
            n_b, 
            block_size, 
            n_repeats, 
            prohibited_pairs.clone()
        ).unwrap();
        
        if result.best_d > min_d {
            min_d = result.best_d;
            best_solution = result;
            best_iter = i;
        }
    }
    Ok((best_solution, best_iter, min_d))
}



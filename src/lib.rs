mod ibd;
mod random_type;


pub mod coincidence_matrix;
pub mod block_array;
pub mod block_result;

use ibd::*;

use block_result::BlockResult;

use anyhow::Error;


/// Finds the best incomplete block design (IBD) by running multiple iterations of the IBD generation algorithm.
///
/// # Arguments
///
/// * `v` - Number of treatments
/// * `n_b` - Number of blocks
/// * `block_size` - Number of treatments per block
/// * `n_repeats` - Number of optimization repeats within each iteration
/// * `iter` - Number of iterations to run
/// * `prohibited_pairs` - Vector of tuples containing pairs of treatments that cannot appear together in any block (0-based indices)
/// * `on_update` - Callback function called after each iteration with the current iteration number and best D-optimality value
///
/// # Returns
///
/// Returns a tuple containing:
/// * The best block design solution found (`BlockResult`)
/// * The iteration number where the best solution was found
/// * The D-optimality value of the best solution
///
/// # Errors
///
/// Returns an error if the IBD generation fails
pub fn find_best_ibd<F>(v: u8, n_b: usize, block_size: usize, n_repeats: usize, iter: usize, prohibited_pairs: Vec<(usize, usize)>, mut on_update: F) -> Result<(BlockResult, usize, f64), Error>
    where
        F: FnMut(usize, f64),
 {
    let mut best_d = 0.0;
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
        
        if result.best_d > best_d {
            best_d = result.best_d;
            best_solution = result;
            best_iter = i;
        }

        on_update(i, best_d);
    }
    Ok((best_solution, best_iter, best_d))
}



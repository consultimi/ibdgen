mod opt_block;

use opt_block::*;
use pretty_print_nalgebra::*;


fn main() {
    let mut min_d = 0.0;
    let mut best_solution = BlockResult::default();
    for _ in 1..10 {
        //let result = opt_block(21,21, 5,5).unwrap();
        let result = opt_block(9,12, 3,5).unwrap();
        //let result = opt_block(7,7, 3,5).unwrap();
        
        //println!("result: {:?}", result);
        if result.best_d > min_d {
            min_d = result.best_d;
            best_solution = result;
        }
    }
    println!("min_d: {:?}", min_d);
    println!("best_solution coincidence: {}", pretty_print!(&best_solution.best_coincidence.coincidence));

    let sorted_block_array = best_solution.best_block_array.as_sorted();
    println!("best_solution block_array: {}", pretty_print!(&sorted_block_array.add_scalar(1)));
    
}

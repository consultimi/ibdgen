mod opt_block;

use nalgebra::DMatrix;
use opt_block::*;
use pretty_print_nalgebra::*;

#[allow(unused)]
fn dm7choose3() -> DMatrix<u8> {
    nalgebra::dmatrix![
        0,0,0,0,0,0;
        1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,1,0,0;
        0,0,0,0,1,0;
        0,0,0,0,0,1
    ]
}

#[allow(unused)]
fn dm9choose3() -> DMatrix<u8> {
    nalgebra::dmatrix![
        0,0,0,0,0,0,0,0;
        1,0,0,0,0,0,0,0;
        0,1,0,0,0,0,0,0;
        0,0,1,0,0,0,0,0;
        0,0,0,1,0,0,0,0;
        0,0,0,0,1,0,0,0;
        0,0,0,0,0,1,0,0;
        0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,1;
    ]
}

#[allow(unused)]
fn dm12choose4() -> DMatrix<u8> {
    nalgebra::dmatrix![
        0,0,0,0,0,0,0,0,0,0,0;
        1,0,0,0,0,0,0,0,0,0,0;
        0,1,0,0,0,0,0,0,0,0,0;
        0,0,1,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0;
        0,0,0,0,1,0,0,0,0,0,0;
        0,0,0,0,0,1,0,0,0,0,0;
        0,0,0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,0,1,0,0;
        0,0,0,0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,0,0,0,1;
    ]
}

fn main() {
    //let x = dm7choose3();
    //let result = opt_block(x.cast::<f64>(), None, 7, vec![3, 3,3,3,3,3,3], 1);
    //assert!(result.is_ok());

    //let x = dm9choose3();
    //let result = opt_block(x.cast::<f64>(), None, 12, vec![3, 3,3,3,3,3,3,3,3,3,3,3], 1);
    //assert!(result.is_ok());

    //let x = dm12choose4();      
    let mut min_d = 0.0;
    let mut best_solution = BlockResult::default();
    for i in 1..10 {
        let result = opt_block(21,21, 5,5).unwrap();
        println!("result: {:?}", result);
        if result.best_d > min_d {
            min_d = result.best_d;
            best_solution = result;
        }
    }
    println!("min_d: {:?}", min_d);
    println!("best_solution coincidence: {}", pretty_print!(&best_solution.best_coincidence.coincidence));
    
}

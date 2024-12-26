mod opt_block;

use nalgebra::DMatrix;
use opt_block::*;

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

    let x = dm12choose4();      
    let result = opt_block(x.cast::<f64>(), None, 33, vec![4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],10 );
    assert!(result.is_ok());
}

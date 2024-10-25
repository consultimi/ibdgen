mod opt_block;

use nalgebra::DMatrix;
use opt_block::*;

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

fn main() {
    let x = dm7choose3();
    let result = opt_block(x.cast::<f64>(), None, 7, vec![3, 3,3,3,3,3,3], true, None, 1, 0);
    assert!(result.is_ok());
}

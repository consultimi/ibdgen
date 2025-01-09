use crate::block_array::BlockArray;
use crate::coincidence_matrix::CoincidenceMatrix;

#[derive(Debug, Default)]
pub struct BlockResult {
    pub best_log_det: f64,
    pub best_block_array: BlockArray,
    pub best_d: f64,
    pub best_diagonality: f64,
    pub best_coincidence: CoincidenceMatrix
}
use crate::block_array::BlockArray;
use crate::coincidence_matrix::CoincidenceMatrix;

/// Represents the result of generating an incomplete block design (IBD).
///
/// This struct holds the best solution found during the IBD generation process,
/// including the block arrangement, various optimality measures, and the coincidence matrix.
///
/// # Fields
///
/// * `best_log_det` - The logarithm of the determinant for the best solution
/// * `best_block_array` - The block arrangement that produced the best result
/// * `best_d` - The D-optimality criterion value for the best solution
/// * `best_diagonality` - A measure of how close the information matrix is to being diagonal
/// * `best_coincidence` - The coincidence matrix for the best solution, showing how often treatments appear together
#[derive(Debug, Default)]
pub struct BlockResult {
    pub best_log_det: f64,
    pub best_block_array: BlockArray,
    pub best_d: f64,
    pub best_diagonality: f64,
    pub best_coincidence: CoincidenceMatrix
}
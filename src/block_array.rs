use nalgebra::DMatrix;

/// Represents a block array, which is a matrix where each row represents a block and each column represents a treatment.
///
/// # Fields
///
/// * `block_array` - The matrix representing the block array
#[derive(Default, Debug, Clone)]
pub struct BlockArray {
    pub block_array: DMatrix<usize>,
}


impl BlockArray {
    pub fn from_block_array(block_array: &DMatrix<usize>) -> Self {
        Self { block_array: block_array.clone() }
    }

    pub fn as_sorted(&mut self) -> DMatrix<usize> {
        let mut block_array_out = self.block_array.clone(); 
        //println!("block_array before sort: {}", pretty_print!(&block_array_out));
        // first, sort each row ascending
        for mut row in block_array_out.row_iter_mut() {
            let mut swapped = true;

            while swapped {
                swapped = false;
                for i in 0..row.len()-1 {
                    if row[i] > row[i + 1] {
                        row.swap_columns(i, i + 1);
                        swapped = true;
                    }
                }
            }
        }

        for _ in 0..block_array_out.ncols() {
            let mut swapped = true;
            while swapped {
                swapped = false;
                for j in 0..block_array_out.nrows()-1 {
                    let row_j = block_array_out.row(j).clone_owned();
                    let row_j1 = block_array_out.row(j+1).clone_owned();
                    let diff = row_j.cast::<i32>() - row_j1.cast::<i32>();
                    //println!("diff: {}", pretty_print!(&diff));
                    if let Some(first_non_zero) = diff.iter().position(|&x| x != 0) {
                        if diff[first_non_zero] > 0 {
                            block_array_out.swap_rows(j, j+1);
                            swapped = true;
                        }
                    }
                }
            }
        }
        // then, sort each column ascending
        block_array_out
    }
}
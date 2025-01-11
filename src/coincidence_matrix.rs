use nalgebra::{DMatrix};

#[derive(Default, Debug, Clone)]
pub struct CoincidenceMatrix {
    pub coincidence: DMatrix<usize>,
}

impl CoincidenceMatrix {


    pub fn r(&self) -> f64 {
        self.coincidence.diagonal().cast::<f64>().mean()
    }

    pub fn variance(&self) -> f64 {
        let mut upper_try_f: nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self.coincidence.upper_triangle().cast::<f64>();
        upper_try_f.fill_lower_triangle_with_upper_triangle();
        upper_try_f.fill_diagonal(self.lambda());
        upper_try_f.variance()
    }

    pub fn lambda(&self) -> f64 {
        let mut upper_try_f: nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self.coincidence.upper_triangle().cast::<f64>();
        upper_try_f.fill_diagonal(0.0);

        let cells = upper_try_f.ncols() * (upper_try_f.ncols() - 1) / 2;
        
        upper_try_f.sum() / (cells as f64)
    }

    pub fn is_bibd(&self) -> bool {
        let diagonal_f = self.coincidence.diagonal().cast::<f64>();
        (diagonal_f.max() == diagonal_f.min()) && (self.variance() == 0.0)
    }

    pub fn from_block_array(block_array: &DMatrix<usize>) -> Self {
        //println!("block_array: {}", pretty_print!(&block_array));   
        //let n = block_array.max() + 1;
        let true_n = block_array.max() + 1;
        let n = block_array.nrows();
        let block_size = block_array.ncols();
        let mut coincidence: DMatrix<usize> = DMatrix::zeros(n, n);
        for block_idx in 0..n {
            let block_elements = block_array.row(block_idx);
            for i in 0..block_size {
                let elem_i = block_elements[i];
                // Diagonal counts total appearances
                coincidence[(elem_i, elem_i)] += 1;
                
                // Upper triangle counts pairwise coincidences
                for j in (i+1)..block_size {
                    
                    let elem_j = block_elements[j];
                    
                    if elem_i < elem_j {
                        coincidence[(elem_i, elem_j)] += 1;
                    } else {
                        coincidence[(elem_j, elem_i)] += 1;
                    }
                }
            }
        }
        Self { coincidence: coincidence.view((0,0),(true_n,true_n)).into() }
    }
}

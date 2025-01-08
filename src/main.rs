mod opt_block;

use opt_block::*;
use pretty_print_nalgebra::*;
use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of varieties
    v: u8,
    
    /// Number of blocks 
    n_b: u8,
    
    /// Size of each block
    block_size: u8,

    /// Number of iterations, that is, how many times the whole process is repeated
    #[arg(long, short, default_value_t = 10)]
    iter: usize,
    
    /// Number of optimization repeats
    #[arg(long, default_value_t = 5)]
    n_repeats: u8,
    
    /// Prohibited pairs as comma-separated values (e.g., "1,2,3,4" means pairs (1,2) and (3,4))
    #[arg(long, short, value_delimiter = ',', num_args = 1..)]
    prohibited_pairs: Option<Vec<String>>,
}

fn main() {
    let args = Args::parse();
    

    let mut prohibited_pairs: Vec<(u8, u8)> = vec![];
    if let Some(pairs) = args.prohibited_pairs {
        if pairs.len() % 2 != 0 {
            println!("Prohibited pairs should be 0 or an even number");
            return;
        }

        // args.prohibited_pairs is a comma-separated list of pairs, each odd and even entry are a pair
        // e.g. "1,2,3,4" -> (1,2) and (3,4)
        for i in (0..pairs.len()).step_by(2) {
            let (a, b) = (pairs[i].parse::<u8>().unwrap(), pairs[i+1].parse::<u8>().unwrap());
            prohibited_pairs.push((a - 1, b - 1));
        }
       

    }
    let mut min_d = 0.0;
    let mut best_solution = BlockResult::default();
    for _ in 1..=args.iter {
        let result = opt_block(
            args.v, 
            args.n_b, 
            args.block_size, 
            args.n_repeats, 
            prohibited_pairs.clone()
        ).unwrap();
        
        if result.best_d > min_d {
            min_d = result.best_d;
            best_solution = result;
        }
    }
    
    println!("Best solution found");
    println!("-------------------");

    let mut coincidence_f = best_solution.best_coincidence.coincidence.clone().cast::<f64>();
    coincidence_f.fill_lower_triangle_with_upper_triangle();
    let udu = coincidence_f.symmetric_eigen();
    //if let Some(udu) = udu {
        println!("eigen: {}", pretty_print!(&udu.eigenvalues));
    //}
    //println!("SINGVAL: {}", pretty_print!(&coincidence_f.singular_values()));
    if best_solution.best_coincidence.is_bibd() {
        println!("Solution is a BIBD with v = {}, k = {}, lambda = {}", args.v, args.block_size, best_solution.best_coincidence.lambda());
    } else {
        println!("Solution is not a BIBD");
    }

    println!("Log Determinant: {:?}", best_solution.best_log_det);
    println!("D-Optimality: {:?}", min_d);
    println!("Diagonality: {:?}", best_solution.best_diagonality);    
    println!("best_solution coincidence: {}", pretty_print!(&best_solution.best_coincidence.coincidence));
    let sorted_block_array = best_solution.best_block_array.as_sorted();
    println!("best_solution block_array: {}", pretty_print!(&sorted_block_array.add_scalar(1)));
    
}

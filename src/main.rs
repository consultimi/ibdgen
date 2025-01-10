use pretty_print_nalgebra::*;
use clap::Parser;


#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of varieties
    v: u8,
    
    /// Number of blocks 
    n_b: usize,
    
    /// Size of each block
    block_size: usize,

    /// Number of iterations, that is, how many times the whole process is repeated
    #[arg(long, short, default_value_t = 10)]
    iter: usize,
    
    /// Number of optimization repeats
    #[arg(long, default_value_t = 5)]
    nrepeats: usize,
    
    /// Prohibited pairs as comma-separated values (e.g., "1,2,3,4" means pairs (1,2) and (3,4))
    #[arg(long, short, value_delimiter = ',', num_args = 1..)]
    prohibited_pairs: Option<Vec<String>>,
}

fn main() {
    let args = Args::parse();
    

    let mut prohibited_pairs: Vec<(usize, usize)> = vec![];
    if let Some(pairs) = args.prohibited_pairs {
        if pairs.len() % 2 != 0 {
            println!("Prohibited pairs should be 0 or an even number");
            return;
        }

        // args.prohibited_pairs is a comma-separated list of pairs, each odd and even entry are a pair
        // e.g. "1,2,3,4" -> (1,2) and (3,4)
        for i in (0..pairs.len()).step_by(2) {
            let (a, b) = (pairs[i].parse::<usize>().unwrap(), pairs[i+1].parse::<usize>().unwrap());
            prohibited_pairs.push((a - 1, b - 1));
        }
    }

    match ibdgen::find_best_ibd(args.v, args.n_b, args.block_size,  args.nrepeats, args.iter, prohibited_pairs) {
        Ok((mut best_solution, best_iter, min_d)) => {
            println!("Best solution found at iteration: {}", best_iter);
            println!("-------------------------------------");
        
            println!("Solution is {} BIBD  (v = {}, k = {}, r = {}, lambda = {}, off-diagonal variance = {})", 
                if best_solution.best_coincidence.is_bibd() { "a" } else { "not a" },
                args.v, args.block_size, best_solution.best_coincidence.r(), best_solution.best_coincidence.lambda(), best_solution.best_coincidence.variance());
        
            println!("Log Determinant: {:?}", best_solution.best_log_det);
            println!("D-Optimality: {:?}", min_d);
            println!("Diagonality: {:?}", best_solution.best_diagonality);    
            println!("best_solution coincidence: {}", pretty_print!(&best_solution.best_coincidence.coincidence));
            let sorted_block_array = best_solution.best_block_array.as_sorted();
            println!("best_solution block_array: {}", pretty_print!(&sorted_block_array.add_scalar(1)));
        
        }
        Err(e) => {
            println!("Error running ibdgen: {:?}", e);
            std::process::exit(1);
        }
    }

}


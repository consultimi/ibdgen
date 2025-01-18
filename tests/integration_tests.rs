use ibdgen::find_best_ibd;

#[test]
fn test_find_best_ibd_basic() {
    let v = 7;
    let n_b = 7; 
    let block_size = 3;
    let n_repeats = 100;
    let iter = 10;
    let prohibited_pairs = vec![];
    
    let result = find_best_ibd(v, n_b, block_size, n_repeats, iter, prohibited_pairs, |_, _| {});
    assert!(result.is_ok());
    
    let (_, best_iter, min_d) = result.unwrap();
    assert!(best_iter <= iter);
    assert!(min_d >= 0.0);
}

#[test]
fn test_find_best_ibd_with_prohibited_pairs() {
    let v = 7;
    let n_b = 7;
    let block_size = 3;
    let n_repeats = 100;
    let iter = 10;
    let prohibited_pairs = vec![(0,1), (2,3)];
    
    let result = find_best_ibd(v, n_b, block_size, n_repeats, iter, prohibited_pairs, |_, _| {});
    assert!(result.is_ok());
    
    let (block_result, _, _) = result.unwrap();
    
    // Check that prohibited pairs don't appear together in any block
    for i in 0..block_result.best_block_array.block_array.nrows() {
        let row = block_result.best_block_array.block_array.row(i);
        let has_01 = row.iter().any(|&x| x == 0) && row.iter().any(|&x| x == 1);
        let has_23 = row.iter().any(|&x| x == 2) && row.iter().any(|&x| x == 3);
        assert!(!(has_01 || has_23));
    }
}

#[test]
fn test_find_best_ibd_with_progress_callback() {
    let v = 7;
    let n_b = 7;
    let block_size = 3;
    let n_repeats = 5;
    let iter = 10;
    let prohibited_pairs = vec![];
    
    let mut iter_update = 0;
    let mut d_update = 0.0;
    let result = find_best_ibd(v, n_b, block_size, n_repeats, iter, prohibited_pairs,  |iter: usize, d: f64| {
        if iter % 2 == 0 {
            iter_update += iter;
            d_update += d;
        }
    });
    assert!(result.is_ok());
    assert_eq!(iter_update, 30);
    assert!(d_update > 0.0);
}

/*
#[test]
fn test_find_best_ibd_invalid_params() {
    let v = 4; // Too small for meaningful BIBD
    let n_b = 3;
    let block_size = 3;
    let n_repeats = 100;
    let iter = 10;
    let prohibited_pairs = vec![];
    
    let result = find_best_ibd(v, n_b, block_size, n_repeats, iter, prohibited_pairs, |_, _| {});
    assert!(result.is_ok()); // Should still complete even with suboptimal parameters
}
 */
const std = @import("std");

const c = @cImport({
    @cInclude("OptBlock.h");
});

fn printMatrix(comptime T: type, array: []const T, rows: usize, cols: usize) void {
    var i: usize = 0;
    while (i < rows) : (i += 1) {
        var j: usize = 0;
        while (j < cols) : (j += 1) {
            std.debug.print("{d:.1} ", .{array[i * cols + j]});
        }
        std.debug.print("\n", .{});
    }
}

pub fn main() !void {

    //double *X,               // Input matrix (N x k)
    //int N,                   // Number of rows
    //int k,                   // Number of columns
    //int nB,                 // Number of blocks
    //int *blocksizes,        // Array of block sizes
    //double *blockFactors,   // Block factors matrix (or NULL)
    //int nRepeats,          // Number of repeats
    //int criterion,         // Optimization criterion
    //bool initRows,         // Whether to use initial rows
    //int *irows            // Initial row assignments (or NULL)

    var X = [_]f64{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    const N = 7;
    const k = 6;
    const nB = 7;
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    const nRepeats = 1;
    const criterion = 0;
    const initRows = false;
    const out = c.BlockOpt(&X, N, k, nB, &blocksizes, null, nRepeats, criterion, initRows, null);
    std.debug.print("out: {}\n", .{out});
    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are s to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    try bw.flush(); // don't forget to flush!
}

test "simple test" {
    var X = [_]f64{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    //printMatrix(f64, &X, 7, 6);
    const N = 7;
    const k = 6;
    const nB = 7;
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    const nRepeats = 1;
    const criterion = 0;
    const initRows = false;
    const out = c.BlockOpt(&X, N, k, nB, &blocksizes, null, nRepeats, criterion, initRows, null);
    try std.testing.expectEqual(@as(i32, 0), out);
}

test "transpose" {
    var X = [_]f64{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    const N = 7;
    const k = 6;
    c.transposeMatrix(&X, N, k);

    const expected = [_]f64{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
    try std.testing.expectEqualSlices(f64, &expected, &X);
}

test "initialize block array" {
    var rows = [_]i32{ 0, 0, 0, 0, 0, 0, 0 };
    var irows = [_]i32{ 1, 2, 3, 1, 2, 3, 1 };
    const N = 7;
    const nB = 7;
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    var block_array = [_]i32{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const expected = [_]i32{ 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    c.initializeBlockArray(&rows, &irows, N, N, nB, &blocksizes, &block_array);
    try std.testing.expectEqualSlices(i32, &block_array, &expected);
}

//void initializeB(
//	int *B,
//	int *rows,
//	int *irows,
//	int N,
//	int Nxb,
//	int nB,
//	int *blocksizes,
//	bool firstRepeat

test "initialize b" {
    var rows = [_]i32{ 0, 0, 0, 0, 0, 0, 0 };
    var irows = [_]i32{ 1, 2, 3, 1, 2, 3, 1 };
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };

    const N = 7;
    const Nxb = 21;
    const nB = 7;
    var B = [_]i32{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const expected = [_]i32{ 0, 2, 4, 6, 3, 1, 5, 0, 4, 3, 5, 6, 2, 1, 0, 3, 6, 1, 5, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    c.initializeB(&B, &rows, &irows, N, Nxb, nB, &blocksizes, true);
    //printMatrix(i32, &B, 7, 3);

    try std.testing.expectEqualSlices(i32, &B, &expected);
}

test "form block means" {
    var X = [_]f64{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
    var B = [_]i32{ 0, 2, 4, 6, 3, 1, 5, 0, 4, 3, 5, 6, 2, 1, 0, 3, 6, 1, 5, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    //const N = 7;
    const k = 6;
    const nB = 7;
    var blockMeans = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    c.formBlockMeans(&blockMeans, &X, &B, k, nB, &blocksizes);
    //printMatrix(f64, &blockMeans, 7, 6);

    const expected = [_]f64{ 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0 };
    try std.testing.expectEqualSlices(f64, &expected, &blockMeans);
}

test "reduce x to t" {
    var X = [_]f64{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
    var B = [_]i32{ 0, 2, 4, 6, 3, 1, 5, 0, 4, 3, 5, 6, 2, 1, 0, 3, 6, 1, 5, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var blockMeans = [_]f64{ 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.0 };
    const k = 6;
    const nB = 7;
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    var blockFactors = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var vec = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var Sc = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var singular = false;
    var T = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    const logDcrit = c.reduceXtoT(&X, &T, &B, &blockMeans, k, nB, &blocksizes, &blockFactors, &vec, &Sc, &singular);
    //printMatrix(f64, &X, 7, 6);
    //printMatrix(i32, &B, 7, 3);
    //printMatrix(f64, &blockMeans, 7, 6);
    //printMatrix(f64, &T, 7, 6);

    try std.testing.expectEqual(@as(f64, 5.4594687633194745e0), logDcrit);
    try std.testing.expectEqual(true, singular);

    const expected = [_]f64{ 2.666666666666667e0, -1.2499999999999999e-1, -1.2499999999999997e-1, -2.4999999999999997e-1, 0, -3.7499999999999994e-1, 3.625e0, -1.1494252873563218e-2, -2.06896551724138e-1, -2.758620689655173e-1, -3.448275862068968e-2, 2.957854406130268e0, -3.6917098445595864e-1, -3.886010362694301e-3, -3.808290155440417e-1, 2.6083765112262522e0, -4.643270981625558e-1, -1.37394471114054e-1, 2.1617281906969046e0, -4.0334890369349363e-1, 1.457445269758617e0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    try std.testing.expectEqualSlices(f64, &expected, &T);
}

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
    var X = [_]f64{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    var B = [_]i32{ 0, 2, 4, 6, 3, 1, 5, 0, 4, 3, 5, 6, 2, 1, 0, 3, 6, 1, 5, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    var blocksizes = [_]i32{ 3, 3, 3, 3, 3, 3, 3 };
    //const N = 7;
    const k = 6;
    const nB = 7;
    var blockMeans = [_]f64{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    c.formBlockMeans(&blockMeans, &X, &B, k, nB, &blocksizes);
    printMatrix(f64, &blockMeans, 7, 6);

    const expected = [_]f64{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
    try std.testing.expectEqualSlices(f64, &blockMeans, &expected);
}

const c = @cImport({
    @cInclude("OptBlock.h");
});
pub fn main() void {
    _ = c.printf("hello\n");
}

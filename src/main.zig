const std = @import("std");
const zibra = @import("zibra.zig");

pub fn main() anyerror!void {
  const Matrix3x3 = zibra.Matrix3x3(f32);
  const m = Matrix3x3.new(&.{ 0, 0, 0, 0, 0, 0, 0, 0, 0 });
  std.debug.print("{}\n", .{ m });
}


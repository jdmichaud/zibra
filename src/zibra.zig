const std = @import("std");
const testing = std.testing;
const fabs = std.math.fabs;
const f32_epsilon = std.math.f32_epsilon;
const Allocator = std.mem.Allocator;

fn typeError(comptime LHST: type, comptime RHST: type, comptime method: []const u8) noreturn {
  switch (@typeInfo(RHST)) {
    .ErrorUnion => @compileError("Cannot " ++ method ++ " " ++ @typeName(LHST) ++ " with an error type. Forgot `try`?"),
    else => @compileError("Cannot " ++ method ++ " " ++ @typeName(LHST) ++ " with " ++ @typeName(RHST)),
  }
}

pub fn Vec(comptime T: type, comptime size: comptime_int) type {
  if (T != f32 and T != f64 and T != f128) {
    @compileError(@typeName(@This()) ++ " expects a float type (f32, f64 or f128) found " ++ @typeName(T));
  }

  if (size != 3 and size != 4) {
    @compileError(@typeName(@This()) ++ " expects a size of 3 or 4, found " ++ @typeName(T));
  }

  return struct {
    const Self = @This();

    data: [size]T,

    pub fn as(self: Self, comptime U: type) Vec(U, size) {
      var newv = Vec(U, size).new(&.{});
      var i: usize = 0;
      while (i < size) : (i += 1) {
        newv.data[i] = @floatCast(U, self.data[i]);
      }
      return newv;
    }

    pub fn new(comptime initializer: []const T) Self {
      const _initializer: [size]T = if (initializer.len == 0)
        [_]T{ 0 } ** size
      else initializer[0..size].*;

      return Self {
        .data = _initializer,
      };
    }

    pub fn copy(self: Self) Self {
      var newv = Vec(T, size).new(&.{});
      var i: usize = 0;
      while (i < size) : (i += 1) {
        newv.data[i] = self.data[i];
      }
      return newv;
    }

    pub fn equal(self: Self, rhs: Self) bool {
      var i: usize = 0;
      while (i < size) : (i += 1) {
        if (fabs(self.data[i] - rhs.data[i]) > f32_epsilon) {
          return false;
        }
      }
      return true;
    }

    pub fn add(self: *Self, rhs: anytype) Self {
      return switch (@TypeOf(rhs)) {
        Self => {
          var i: usize = 0;
          while (i < size) : (i += 1) {
            self.data[i] += rhs.data[i];
          }
          return self.*;
        },
        comptime_int, comptime_float => {
          var i: usize = 0;
          while (i < size) : (i += 1) {
            self.data[i] += rhs;
          }
          return self.*;
        },
        else => |t| typeError(T, t, "add"),
      };
    }

    pub fn sub(self: *Self, rhs: anytype) Self {
      return switch (@TypeOf(rhs)) {
        Self => {
          var i: usize = 0;
          while (i < size) : (i += 1) {
            self.data[i] -= rhs.data[i];
          }
          return self.*;
        },
        comptime_int, comptime_float => {
          var i: usize = 0;
          while (i < size) : (i += 1) {
            self.data[i] -= rhs;
          }
          return self.*;
        },
        else => |t| typeError(T, t, "add"),
      };
    }

    pub fn dot(self: Self, rhs: Self) T {
      var result: T = 0;
      var i: usize = 0;
      while (i < size) : (i += 1) {
        result += self.data[i] * rhs.data[i];
      }
      return result;
    }

    pub fn cross(self: *Self, rhs: Self) Self {
      if (size != 3) {
        @compileError(@typeName(@This()) ++ ": cross produce implemented only for vector of size " ** size);
      }
      var x = self.data[1] * rhs.data[2] - self.data[2] * rhs.data[1];
      var y = self.data[2] * rhs.data[0] - self.data[0] * rhs.data[2];
      var z = self.data[0] * rhs.data[1] - self.data[1] * rhs.data[0];
      self.data[0] = x; self.data[1] = y; self.data[2] = z;
      return self.*;
    }
  };

}

test "Vec" {
  const V3 = Vec(f64, 3);
  var v = V3.new(&.{});
  try testing.expect(v.equal(V3.new(&.{ 0, 0, 0 })));
  _ = v.add(V3.new(&.{ 1, 2, 3 }));
  try testing.expect(v.equal(V3.new(&.{ 1, 2, 3 })));
  _ = v.add(2.0);
  try testing.expect(v.equal(V3.new(&.{ 3, 4, 5 })));
  try testing.expect(v.copy().sub(1).equal(V3.new(&.{ 2, 3, 4 })));
  try testing.expect(v.equal(V3.new(&.{ 3, 4, 5 })));
  try testing.expect(V3.new(&.{ 1, 2, 3 }).dot(V3.new(&.{ 1, 2, 3 })) == 14);
  try testing.expect(V3.new(&.{ 1, 2, 3 }).as(f32).equal(Vec(f32, 3).new(&.{ 1, 2, 3 })));
  try testing.expect(V3.new(&.{ 1, 0, 0 }).cross(V3.new(&.{ 0, 1, 0 })).equal(V3.new(&.{ 0, 0, 1 })));
  try testing.expect(V3.new(&.{ 0, 1, 0 }).cross(V3.new(&.{ 0, 0, 1 })).equal(V3.new(&.{ 1, 0, 0 })));
  try testing.expect(V3.new(&.{ 0, 0, 1 }).cross(V3.new(&.{ 1, 0, 0 })).equal(V3.new(&.{ 0, 1, 0 })));
}

pub fn Vec3(comptime T: type) type {
  if (T != f32 and T != f64 and T != f128) {
    @compileError(@typeName(@This()) ++ " expects a float type (f32, f64 or f128) found " ++ @typeName(T));
  }

  return struct {
    const Self = @This();

    x: T,
    y: T,
    z: T,

    pub fn as(self: Self, comptime U: type) Vec3(U) {
      var newv = Vec3(U).new(&.{});
      newv.x = @floatCast(U, self.x);
      newv.y = @floatCast(U, self.y);
      newv.z = @floatCast(U, self.z);
      return newv;
    }

    pub fn new(comptime initializer: []const T) Self {
      const _initializer = if (initializer.len == 0)
        &[_]T{ 0, 0, 0 }
      else initializer;

      return Self {
        .x = _initializer[0],
        .y = _initializer[1],
        .z = _initializer[2],
      };
    }

    pub fn copy(self: Self) Self {
      var newv = Self.new(&.{});
      newv.x = self.x;
      newv.y = self.y;
      newv.z = self.z;
      return newv;
    }

    pub fn equal(self: Self, rhs: Self) bool {
      return fabs(self.x - rhs.x) < f32_epsilon and
             fabs(self.y - rhs.y) < f32_epsilon and
             fabs(self.z - rhs.z) < f32_epsilon;
    }

    pub fn add(self: *Self, rhs: anytype) Self {
      return switch (@TypeOf(rhs)) {
        Self => {
          self.x += rhs.x;
          self.y += rhs.y;
          self.z += rhs.z;
          return self.*;
        },
        comptime_int, comptime_float => {
          self.x += rhs;
          self.y += rhs;
          self.z += rhs;
          return self.*;
        },
        else => |t| typeError(T, t, "add"),
      };
    }

    pub fn sub(self: *Self, rhs: anytype) Self {
      return switch (@TypeOf(rhs)) {
        Self => {
          self.x -= rhs.x;
          self.y -= rhs.y;
          self.z -= rhs.z;
          return self.*;
        },
        comptime_int, comptime_float => {
          self.x -= rhs;
          self.y -= rhs;
          self.z -= rhs;
          return self.*;
        },
        else => |t| typeError(T, t, "add"),
      };
    }

    pub fn dot(self: Self, rhs: Self) T {
      return self.x * rhs.x +
             self.y * rhs.y +
             self.z * rhs.z;
    }

    pub fn cross(self: *Self, rhs: Self) Self {
      var x = self.y * rhs.z - self.z * rhs.y;
      var y = self.z * rhs.x - self.x * rhs.z;
      var z = self.x * rhs.y - self.y * rhs.x;
      self.x = x; self.y = y; self.z = z;
      return self.*;
    }
  };
}

test "Vec3" {
  const V3 = Vec3(f64);
  var v = V3.new(&.{});
  try testing.expect(v.equal(V3.new(&.{ 0, 0, 0 })));
  _ = v.add(V3.new(&.{ 1, 2, 3 }));
  try testing.expect(v.equal(V3.new(&.{ 1, 2, 3 })));
  _ = v.add(2.0);
  try testing.expect(v.equal(V3.new(&.{ 3, 4, 5 })));
  try testing.expect(v.copy().sub(1).equal(V3.new(&.{ 2, 3, 4 })));
  try testing.expect(v.equal(V3.new(&.{ 3, 4, 5 })));
  try testing.expect(V3.new(&.{ 1, 2, 3 }).dot(V3.new(&.{ 1, 2, 3 })) == 14);
  try testing.expect(V3.new(&.{ 1, 2, 3 }).as(f32).equal(Vec3(f32).new(&.{ 1, 2, 3 })));
  try testing.expect(V3.new(&.{ 1, 0, 0 }).cross(V3.new(&.{ 0, 1, 0 })).equal(V3.new(&.{ 0, 0, 1 })));
  try testing.expect(V3.new(&.{ 0, 1, 0 }).cross(V3.new(&.{ 0, 0, 1 })).equal(V3.new(&.{ 1, 0, 0 })));
  try testing.expect(V3.new(&.{ 0, 0, 1 }).cross(V3.new(&.{ 1, 0, 0 })).equal(V3.new(&.{ 0, 1, 0 })));
}

pub fn Matrix(comptime T: type, comptime size: comptime_int) type {
  if (T != f32 and T != f64 and T != f128) {
    @compileError(@typeName(@This()) ++ " expects a float type (f32, f64 or f128) found " ++ @typeName(T));
  }

  if (size != 3 and size != 4) {
    @compileError(@typeName(@This()) ++ " expects a size of 3 or 4, found " ++ @typeName(T));
  }

  return struct {
    const Self = @This();

    data: [size * size]T,

    pub fn as(self: Self, comptime U: type) Matrix(U, size) {
      var newv = Matrix(U, size).new(&.{});
      var i: usize = 0;
      while (i < size) : (i += 1) {
        newv.data[i] = @floatCast(U, self.data[i]);
      }
      return newv;
    }

    pub fn new(comptime initializer: []const T) Self {
      const _initializer: [size]T = if (initializer.len == 0)
        [_]T{ 0 } ** size
      else initializer[0..size].*;

      return Self {
        .data = _initializer,
      };
    }

    pub fn copy(self: Self) Self {
      var newv = Matrix(T, size).new(&.{});
      var i: usize = 0;
      while (i < size) : (i += 1) {
        newv.data[i] = self.data[i];
      }
      return newv;
    }

    pub fn equal(self: Self, rhs: Self) bool {
      var i: usize = 0;
      while (i < size) : (i += 1) {
        if (fabs(self.data[i] - rhs.data[i]) > f32_epsilon) {
          return false;
        }
      }
      return true;
    }

    pub fn identity() Self {
      return if (size == 3)
        @This().new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })
      else
        @This().new(&.{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 });
    }

    pub fn add(self: *Self, rhs: Self) Self {
      for (self.data) |_, i| {
        self.data[i] += rhs.data[i];
      }
      return self.*;
    }

    pub fn neg(self: Self) Self {
      for (self.data) |_, i| {
        self.data[i] -= self.data[i];
      }
      return self;
    }

    pub fn sub(self: *Self, rhs: Self) Self {
      for (self.data) |_, i| {
        self.data[i] -= rhs.data[i];
      }
      return self.*;
    }

    pub fn mul(self: *Self, rhs: anytype) @TypeOf(rhs) {
      switch (@TypeOf(rhs)) {
        Self => {
          var tmp: [size * size]T = [_]T{0} ** (size * size);
          tmp[0] = self.data[0] * rhs.data[0] + self.data[1] * rhs.data[3] + self.data[2] * rhs.data[6];
          tmp[1] = self.data[0] * rhs.data[1] + self.data[1] * rhs.data[4] + self.data[2] * rhs.data[7];
          tmp[2] = self.data[0] * rhs.data[2] + self.data[1] * rhs.data[5] + self.data[2] * rhs.data[8];
          tmp[3] = self.data[3] * rhs.data[0] + self.data[4] * rhs.data[3] + self.data[5] * rhs.data[6];
          tmp[4] = self.data[3] * rhs.data[1] + self.data[4] * rhs.data[4] + self.data[5] * rhs.data[7];
          tmp[5] = self.data[3] * rhs.data[2] + self.data[4] * rhs.data[5] + self.data[5] * rhs.data[8];
          tmp[6] = self.data[6] * rhs.data[0] + self.data[7] * rhs.data[3] + self.data[8] * rhs.data[6];
          tmp[7] = self.data[6] * rhs.data[1] + self.data[7] * rhs.data[4] + self.data[8] * rhs.data[7];
          tmp[8] = self.data[6] * rhs.data[2] + self.data[7] * rhs.data[5] + self.data[8] * rhs.data[8];
          std.mem.copy(T, self.data[0..], tmp[0..]);
          return self.*;
        },
        Vec3(T) => {
          var result = Vec3(T).new(&.{});
          result.x = self.data[0] * rhs.x + self.data[1] * rhs.y + self.data[2] * rhs.z;
          result.y = self.data[3] * rhs.x + self.data[4] * rhs.y + self.data[5] * rhs.z;
          result.z = self.data[6] * rhs.x + self.data[7] * rhs.y + self.data[8] * rhs.z;
          return result;
        },
        else => |t| typeError(T, t, "mul"),
      }
    }

    //     a b c
    // det(d e f) = a * (ei - hf) - b * (di - gf) + c * (dh - ge)
    //     g h i
    fn determinant3(self: Self) T {
      return self.data[0] * (self.data[4] * self.data[8] - self.data[7] * self.data[5])
           - self.data[1] * (self.data[3] * self.data[8] - self.data[6] * self.data[5])
           + self.data[2] * (self.data[3] * self.data[7] - self.data[6] * self.data[4])
      ;
    }

    //
    //     a b c d
    //     e f g h
    // det(i j k l) =   a * (f * (kp - ol) - g * (jp - nl) + h * (jo - nk))
    //     m n o p    - b * (e * (kp - ol) - g * (ip - ml) + h * (io - mk))
    //                + c * (e * (jp - nl) - f * (ip - ml) + h * (in - mj))
    //                - d * (e * (jo - nk) - f * (io - mk) + g * (in - mj))
    //
    fn determinant4(self: Self) T {
      return
          self.data[0] * (self.data[5] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[6] * (self.data[9] * self.data[15] - self.data[13] * self.data[11]) + self.data[7] * (self.data[9] * self.data[14] - self.data[13] * self.data[10]))
        - self.data[1] * (self.data[4] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[6] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[7] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]))
        + self.data[2] * (self.data[4] * ( self.data[9] * self.data[15] - self.data[13] * self.data[11]) - self.data[5] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[7] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]))
        - self.data[3] * (self.data[4] * ( self.data[9] * self.data[14] - self.data[13] * self.data[10]) - self.data[5] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]) + self.data[6] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]))
        ;
    }

    pub fn determinant(self: Self) T {
      if (size == 3) {
        return self.determinant3();
      } else {
        if (size == 4) {
          return self.determinant4();
        }
      }
      @compileError(@typeName(@This()) ++ ": determinant not implemented for matrix of size " ++ size);
    }

    // To compute an inverse:
    //
    // 1. First check the determinant is not zero
    // 2. Compute the cofactor matrix
    // 3. Transpose the cofactor matrix to get the adjugate matrix
    // 4. Divide the adjugate matrix by the determinant
    pub fn inverse3(self: *Self) !Self {
      // Compute determinant to check if matrix is invertible
      const det: T = self.determinant();
      if (fabs(det) < f32_epsilon) return error.NoInverse;
      // Compute adjugate matrix (transpose of the cofactors)
      var adjugate = [_]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      adjugate[0] =  (self.data[4] * self.data[8] - self.data[7] * self.data[5]);
      adjugate[3] = -(self.data[3] * self.data[8] - self.data[6] * self.data[5]);
      adjugate[6] =  (self.data[3] * self.data[7] - self.data[6] * self.data[4]);
      adjugate[1] = -(self.data[1] * self.data[8] - self.data[7] * self.data[2]);
      adjugate[4] =  (self.data[0] * self.data[8] - self.data[6] * self.data[2]);
      adjugate[7] = -(self.data[0] * self.data[7] - self.data[6] * self.data[1]);
      adjugate[2] =  (self.data[1] * self.data[5] - self.data[4] * self.data[2]);
      adjugate[5] = -(self.data[0] * self.data[5] - self.data[3] * self.data[2]);
      adjugate[8] =  (self.data[0] * self.data[4] - self.data[3] * self.data[1]);
      // Divide the cofactor matrix by the dererminant
      const idet: T = 1 / det;
      self.data[0] = adjugate[0] * idet;
      self.data[1] = adjugate[1] * idet;
      self.data[2] = adjugate[2] * idet;
      self.data[3] = adjugate[3] * idet;
      self.data[4] = adjugate[4] * idet;
      self.data[5] = adjugate[5] * idet;
      self.data[6] = adjugate[6] * idet;
      self.data[7] = adjugate[7] * idet;
      self.data[8] = adjugate[8] * idet;

      return self.*;
    }

    pub fn inverse4(self: *Self) !Self {
      // Compute determinant to check if matrix is invertible
      const det = determinant4(m);
      if (fabs(det) < f32_epsilon) return error.NoInverse;
      // Compute adjugate matrix (transpose of the cofactors)
      //          a b c d     A11 A12 A13 A14
      // cofactor(e f g h) = (A21 A22 A23 A24)
      //          i j k l     A31 A32 A33 A34
      //          m n o p     A41 A42 A43 A44
      //
      //   A11 =  a * (f * (k * p - o * l) - g * (j * p - n * l) + h * (j * o - n * k))
      //   A12 = -b * (e * (k * p - o * l) - g * (i * p - m * l) + h * (i * o - m * k))
      //   A13 =  c * (e * (j * p - n * l) - f * (i * p - m * l) + h * (i * n - m * j))
      //   A14 = -d * (e * (j * o - n * k) - f * (i * o - m * k) + g * (i * n - m * j))
      //   A21 =  e * (b * (k * p - o * l) - c * (j * p - n * l) + d * (j * o - n * k))
      //   A22 = -f * (a * (k * p - o * l) - c * (i * p - m * l) + d * (i * o - m * k))
      //   A23 =  g * (a * (j * p - n * l) - b * (i * p - m * l) + d * (i * n - m * j))
      //   A24 = -h * (a * (j * o - n * k) - b * (i * o - m * k) + c * (i * n - m * j))
      //   A31 =  i * (b * (g * p - o * h) - c * (f * p - n * h) + d * (f * o - n * g))
      //   A32 = -j * (a * (g * p - o * h) - c * (e * p - m * h) + d * (e * o - m * g))
      //   A33 =  k * (a * (f * p - n * h) - b * (e * p - m * h) + d * (e * n - m * f))
      //   A34 = -l * (a * (f * o - n * g) - b * (e * o - m * g) + c * (e * n - m * f))
      //   A41 =  m * (b * (g * l - k * h) - c * (f * l - j * h) + d * (f * k - j * g))
      //   A42 = -n * (a * (g * l - k * h) - c * (e * l - i * h) + d * (e * k - i * g))
      //   A43 =  o * (a * (f * l - j * h) - b * (e * l - i * h) + d * (e * j - i * f))
      //   A44 = -p * (a * (f * k - g * j) - b * (e * k - i * g) + c * (e * j - i * f))
      var adjugate = [_]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      adjugate[0]  =  (self.data[5] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[6] * (self.data[9] * self.data[15] - self.data[13] * self.data[11]) + self.data[7] * (self.data[9] * self.data[14] - self.data[13] * self.data[10]));
      adjugate[4]  = -(self.data[4] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[6] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[7] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]));
      adjugate[8]  =  (self.data[4] * ( self.data[9] * self.data[15] - self.data[13] * self.data[11]) - self.data[5] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[7] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]));
      adjugate[12] = -(self.data[4] * ( self.data[9] * self.data[14] - self.data[13] * self.data[10]) - self.data[5] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]) + self.data[6] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]));
      adjugate[1]  = -(self.data[1] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[2] * (self.data[9] * self.data[15] - self.data[13] * self.data[11]) + self.data[3] * (self.data[9] * self.data[14] - self.data[13] * self.data[10]));
      adjugate[5]  =  (self.data[0] * (self.data[10] * self.data[15] - self.data[14] * self.data[11]) - self.data[2] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[3] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]));
      adjugate[9]  = -(self.data[0] * ( self.data[9] * self.data[15] - self.data[13] * self.data[11]) - self.data[1] * (self.data[8] * self.data[15] - self.data[12] * self.data[11]) + self.data[3] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]));
      adjugate[13] =  (self.data[0] * ( self.data[9] * self.data[14] - self.data[13] * self.data[10]) - self.data[1] * (self.data[8] * self.data[14] - self.data[12] * self.data[10]) + self.data[2] * (self.data[8] * self.data[13] - self.data[12] *  self.data[9]));
      adjugate[2]  =  (self.data[1] * ( self.data[6] * self.data[15] - self.data[14] *  self.data[7]) - self.data[2] * (self.data[5] * self.data[15] - self.data[13] *  self.data[7]) + self.data[3] * (self.data[5] * self.data[14] - self.data[13] *  self.data[6]));
      adjugate[6]  = -(self.data[0] * ( self.data[6] * self.data[15] - self.data[14] *  self.data[7]) - self.data[2] * (self.data[4] * self.data[15] - self.data[12] *  self.data[7]) + self.data[3] * (self.data[4] * self.data[14] - self.data[12] *  self.data[6]));
      adjugate[10] =  (self.data[0] * ( self.data[5] * self.data[15] - self.data[13] *  self.data[7]) - self.data[1] * (self.data[4] * self.data[15] - self.data[12] *  self.data[7]) + self.data[3] * (self.data[4] * self.data[13] - self.data[12] *  self.data[5]));
      adjugate[14] = -(self.data[0] * ( self.data[5] * self.data[14] - self.data[13] *  self.data[6]) - self.data[1] * (self.data[4] * self.data[14] - self.data[12] *  self.data[6]) + self.data[2] * (self.data[4] * self.data[13] - self.data[12] *  self.data[5]));
      adjugate[3]  = -(self.data[1] * ( self.data[6] * self.data[11] - self.data[10] *  self.data[7]) - self.data[2] * (self.data[5] * self.data[11] -  self.data[9] *  self.data[7]) + self.data[3] * (self.data[5] * self.data[10] -  self.data[9] *  self.data[6]));
      adjugate[7]  =  (self.data[0] * ( self.data[6] * self.data[11] - self.data[10] *  self.data[7]) - self.data[2] * (self.data[4] * self.data[11] -  self.data[8] *  self.data[7]) + self.data[3] * (self.data[4] * self.data[10] -  self.data[8] *  self.data[6]));
      adjugate[11] = -(self.data[0] * ( self.data[5] * self.data[11] -  self.data[9] *  self.data[7]) - self.data[1] * (self.data[4] * self.data[11] -  self.data[8] *  self.data[7]) + self.data[3] * (self.data[4] *  self.data[9] -  self.data[8] *  self.data[5]));
      adjugate[15] =  (self.data[0] * ( self.data[5] * self.data[10] -  self.data[6] *  self.data[9]) - self.data[1] * (self.data[4] * self.data[10] -  self.data[8] *  self.data[6]) + self.data[2] * (self.data[4] *  self.data[9] -  self.data[8] *  self.data[5]));
      // Divide the cofactor matrix by the dererminant
      const idet = 1 / det;
      self.data[0] = adjugate[0] * idet;
      self.data[1] = adjugate[1] * idet;
      self.data[2] = adjugate[2] * idet;
      self.data[3] = adjugate[3] * idet;
      self.data[4] = adjugate[4] * idet;
      self.data[5] = adjugate[5] * idet;
      self.data[6] = adjugate[6] * idet;
      self.data[7] = adjugate[7] * idet;
      self.data[8] = adjugate[8] * idet;
      self.data[9] = adjugate[9] * idet;
      self.data[10] = adjugate[10] * idet;
      self.data[11] = adjugate[11] * idet;
      self.data[12] = adjugate[12] * idet;
      self.data[13] = adjugate[13] * idet;
      self.data[14] = adjugate[14] * idet;
      self.data[15] = adjugate[15] * idet;

      return self.*;
    }

    pub fn inverse(self: *Self) !Self {
      if (size == 3) {
        return self.inverse3();
      } else {
        if (size == 4) {
          return self.inverse4();
        }
      }
      @compileError(@typeName(@This()) ++ ": inverse not implemented for matrix of size " ++ size);
    }

  };
}

test "Matrix" {
  const M3 = Matrix(f32, 3);
  var m = M3.new(&.{});
  try testing.expect(m.equal(M3.new(&.{ 0, 0, 0, 0, 0, 0, 0, 0, 0 })));
  var id = M3.identity();
  try testing.expect(!m.equal(id));
  try testing.expect(id.equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  try testing.expect(id.copy().equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  try testing.expect(id.copy().add(id).equal(M3.new(&.{ 2, 0, 0, 0, 2, 0, 0, 0, 2 })));
  try testing.expect(id.equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  _ = m.add(id);
  try testing.expect(m.equal(id));
  _ = m.sub(id);
  try testing.expect(m.equal(M3.new(&.{ 0, 0, 0, 0, 0, 0, 0, 0, 0 })));
  try testing.expectEqual(id.determinant(), 1);
  try testing.expect((try id.inverse()).equal(id));
  // try testing.expectError(error.NoInverse, M3.new(&.{}).inverse());
  try testing.expect(
    M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .mul(M3.new(&.{ 9, 8, 7, 6, 5, 4, 3, 2, 1 }))
    .equal(M3.new(&.{ 30, 24, 18, 84, 69, 54, 138, 114, 90 })));
  const V3 = Vec3(f32);
  try testing.expect(
    M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .mul(V3.new(&.{ 3, 2, 1 }))
    .equal(V3.new(&.{ 10, 28, 46 })));
  try testing.expect(M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .as(f128)
    .equal(Matrix(f128, 3).new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })));
}

pub fn Matrix3x3(comptime T: type) type {
  if (T != f32 and T != f64 and T != f128) {
    @compileError(@typeName(@This()) ++ " expects a float type (f32, f64 or f128) found " ++ @typeName(T));
  }

  return struct {
    const Self = @This();

    data: [9]T,

    pub fn as(self: Self, comptime U: type) Matrix3x3(U) {
      var newm = Matrix3x3(U).new(&.{});
      for (newm.data) |_, i| {
        newm.data[i] = @floatCast(U, self.data[i]);
      }
      return newm;
    }

    pub fn new(comptime initializer: []const T) Self {
      var _initializer: [9]T = if (initializer.len == 0)
        [_]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0 }
      else initializer[0..9].*;

      return Self {
        .data = _initializer,
      };
    }

    pub fn copy(self: Self) Self {
      var newm = Self.new(&.{});
      std.mem.copy(T, newm.data[0..], self.data[0..]);
      return newm;
    }

    pub fn identity() Self {
      return @This().new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 });
    }

    //     a b c
    // det(d e f) = a * (ei - hf) - b * (di - gf) + c * (dh - ge)
    //     g h i
    pub fn determinant(self: Self) T {
      return self.data[0] * (self.data[4] * self.data[8] - self.data[7] * self.data[5])
           - self.data[1] * (self.data[3] * self.data[8] - self.data[6] * self.data[5])
           + self.data[2] * (self.data[3] * self.data[7] - self.data[6] * self.data[4])
      ;
    }

    // To compute an inverse:
    //
    // 1. First check the determinant is not zero
    // 2. Compute the cofactor matrix
    // 3. Transpose the cofactor matrix to get the adjugate matrix
    // 4. Divide the adjugate matrix by the determinant
    pub fn inverse(self: *Self) !Self {
      // Compute determinant to check if matrix is invertible
      const det: T = self.determinant();
      if (fabs(det) < f32_epsilon) return error.NoInverse;
      // Compute adjugate matrix (transpose of the cofactors)
      var adjugate = [_]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      adjugate[0] =  (self.data[4] * self.data[8] - self.data[7] * self.data[5]);
      adjugate[3] = -(self.data[3] * self.data[8] - self.data[6] * self.data[5]);
      adjugate[6] =  (self.data[3] * self.data[7] - self.data[6] * self.data[4]);
      adjugate[1] = -(self.data[1] * self.data[8] - self.data[7] * self.data[2]);
      adjugate[4] =  (self.data[0] * self.data[8] - self.data[6] * self.data[2]);
      adjugate[7] = -(self.data[0] * self.data[7] - self.data[6] * self.data[1]);
      adjugate[2] =  (self.data[1] * self.data[5] - self.data[4] * self.data[2]);
      adjugate[5] = -(self.data[0] * self.data[5] - self.data[3] * self.data[2]);
      adjugate[8] =  (self.data[0] * self.data[4] - self.data[3] * self.data[1]);
      // Divide the cofactor matrix by the dererminant
      const idet: T = 1 / det;
      self.data[0] = adjugate[0] * idet;
      self.data[1] = adjugate[1] * idet;
      self.data[2] = adjugate[2] * idet;
      self.data[3] = adjugate[3] * idet;
      self.data[4] = adjugate[4] * idet;
      self.data[5] = adjugate[5] * idet;
      self.data[6] = adjugate[6] * idet;
      self.data[7] = adjugate[7] * idet;
      self.data[8] = adjugate[8] * idet;

      return self.*;
    }

    pub fn equal(self: Self, rhs: Self) bool {
      for (self.data) |_, i| {
        if (fabs(self.data[i] - rhs.data[i]) > f32_epsilon)
          return false;
      }
      return true;
    }

    pub fn add(self: *Self, rhs: Self) Self {
      for (self.data) |_, i| {
        self.data[i] += rhs.data[i];
      }
      return self.*;
    }

    pub fn neg(self: Self) Self {
      for (self.data) |_, i| {
        self.data[i] -= self.data[i];
      }
      return self;
    }

    pub fn sub(self: *Self, rhs: Self) Self {
      for (self.data) |_, i| {
        self.data[i] -= rhs.data[i];
      }
      return self.*;
    }

    pub fn mul(self: *Self, rhs: anytype) @TypeOf(rhs) {
      switch (@TypeOf(rhs)) {
        Self => {
          var tmp: [9]T = [_]T{0} ** 9;
          tmp[0] = self.data[0] * rhs.data[0] + self.data[1] * rhs.data[3] + self.data[2] * rhs.data[6];
          tmp[1] = self.data[0] * rhs.data[1] + self.data[1] * rhs.data[4] + self.data[2] * rhs.data[7];
          tmp[2] = self.data[0] * rhs.data[2] + self.data[1] * rhs.data[5] + self.data[2] * rhs.data[8];
          tmp[3] = self.data[3] * rhs.data[0] + self.data[4] * rhs.data[3] + self.data[5] * rhs.data[6];
          tmp[4] = self.data[3] * rhs.data[1] + self.data[4] * rhs.data[4] + self.data[5] * rhs.data[7];
          tmp[5] = self.data[3] * rhs.data[2] + self.data[4] * rhs.data[5] + self.data[5] * rhs.data[8];
          tmp[6] = self.data[6] * rhs.data[0] + self.data[7] * rhs.data[3] + self.data[8] * rhs.data[6];
          tmp[7] = self.data[6] * rhs.data[1] + self.data[7] * rhs.data[4] + self.data[8] * rhs.data[7];
          tmp[8] = self.data[6] * rhs.data[2] + self.data[7] * rhs.data[5] + self.data[8] * rhs.data[8];
          std.mem.copy(T, self.data[0..], tmp[0..]);
          return self.*;
        },
        Vec3(T) => {
          var result = Vec3(T).new(&.{});
          result.x = self.data[0] * rhs.x + self.data[1] * rhs.y + self.data[2] * rhs.z;
          result.y = self.data[3] * rhs.x + self.data[4] * rhs.y + self.data[5] * rhs.z;
          result.z = self.data[6] * rhs.x + self.data[7] * rhs.y + self.data[8] * rhs.z;
          return result;
        },
        else => |t| typeError(T, t, "mul"),
      }
    }
  };
}

test "Matrix3x3" {
  const M3 = Matrix3x3(f32);
  var m = M3.new(&.{});
  try testing.expect(m.equal(M3.new(&.{ 0, 0, 0, 0, 0, 0, 0, 0, 0 })));
  var id = M3.identity();
  try testing.expect(!m.equal(id));
  try testing.expect(id.equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  try testing.expect(id.copy().equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  try testing.expect(id.copy().add(id).equal(M3.new(&.{ 2, 0, 0, 0, 2, 0, 0, 0, 2 })));
  try testing.expect(id.equal(M3.new(&.{ 1, 0, 0, 0, 1, 0, 0, 0, 1 })));
  _ = m.add(id);
  try testing.expect(m.equal(id));
  _ = m.sub(id);
  try testing.expect(m.equal(M3.new(&.{ 0, 0, 0, 0, 0, 0, 0, 0, 0 })));
  try testing.expectEqual(id.determinant(), 1);
  try testing.expect((try id.inverse()).equal(id));
  // try testing.expectError(error.NoInverse, M3.new(&.{}).inverse());
  try testing.expect(
    M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .mul(M3.new(&.{ 9, 8, 7, 6, 5, 4, 3, 2, 1 }))
    .equal(M3.new(&.{ 30, 24, 18, 84, 69, 54, 138, 114, 90 })));
  const V3 = Vec3(f32);
  try testing.expect(
    M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .mul(V3.new(&.{ 3, 2, 1 }))
    .equal(V3.new(&.{ 10, 28, 46 })));
  try testing.expect(M3.new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })
    .as(f128)
    .equal(Matrix3x3(f128).new(&.{ 1, 2, 3, 4, 5, 6, 7, 8, 9 })));
}


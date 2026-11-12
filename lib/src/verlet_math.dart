import 'dart:math' as math;

/// Immutable 3D vector for XPBD soft-body & Verlet physics.
class VerletVec3 {
  const VerletVec3(this.x, this.y, this.z);

  final double x;
  final double y;
  final double z;

  static const VerletVec3 zero = VerletVec3(0, 0, 0);

  VerletVec3 operator +(VerletVec3 other) =>
      VerletVec3(x + other.x, y + other.y, z + other.z);

  VerletVec3 operator -(VerletVec3 other) =>
      VerletVec3(x - other.x, y - other.y, z - other.z);

  VerletVec3 operator *(double s) => VerletVec3(x * s, y * s, z * s);

  VerletVec3 operator -() => VerletVec3(-x, -y, -z);

  double dot(VerletVec3 other) => x * other.x + y * other.y + z * other.z;

  VerletVec3 cross(VerletVec3 other) => VerletVec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x,
      );

  double get lengthSquared => x * x + y * y + z * z;

  double get length => math.sqrt(lengthSquared);

  VerletVec3 get normalized {
    final len = length;
    if (len <= 1e-12) return zero;
    final inv = 1.0 / len;
    return VerletVec3(x * inv, y * inv, z * inv);
  }

  double distanceTo(VerletVec3 other) => (this - other).length;

  bool get isFinite => x.isFinite && y.isFinite && z.isFinite;

  @override
  bool operator ==(Object other) =>
      identical(this, other) ||
      other is VerletVec3 && x == other.x && y == other.y && z == other.z;

  @override
  int get hashCode => Object.hash(x, y, z);

  @override
  String toString() => 'VerletVec3($x, $y, $z)';
}

/// Mutable particle state for XPBD / Verlet integration.
class VerletNode {
  VerletNode({
    required this.position,
    this.invMass = 1.0,
    VerletVec3? velocity,
    VerletVec3? prevPosition,
  })  : velocity = velocity ?? VerletVec3.zero,
        prevPosition = prevPosition ?? position;

  VerletVec3 position;
  VerletVec3 prevPosition;
  VerletVec3 velocity;
  double invMass;

  bool get isPinned => invMass <= 0.0;
}

import 'verlet_math.dart';

double _solveDistancePair(
  List<VerletNode> nodes,
  int indexA,
  int indexB,
  double restLength,
  double compliance,
  double lambda,
  double invDtSq,
) {
  final a = nodes[indexA];
  final b = nodes[indexB];
  final wA = a.isPinned ? 0.0 : a.invMass;
  final wB = b.isPinned ? 0.0 : b.invMass;
  final wSum = wA + wB;
  if (wSum <= 0.0) return lambda;

  final diff = b.position - a.position;
  final dist = diff.length;
  if (dist <= 1e-12) return lambda;

  final c = dist - restLength;
  final alphaTilde = compliance * invDtSq;
  final denom = wSum + alphaTilde;
  if (denom <= 1e-24) return lambda;

  final dLambda = (-c - alphaTilde * lambda) / denom;
  final nextLambda = lambda + dLambda;

  final n = diff * (1.0 / dist);
  if (wA > 0.0) {
    a.position = a.position - n * (wA * dLambda);
  }
  if (wB > 0.0) {
    b.position = b.position + n * (wB * dLambda);
  }
  return nextLambda;
}

/// XPBD distance constraint between two particle indices.
class VerletDistanceLink {
  VerletDistanceLink(
    this.nodeA,
    this.nodeB,
    this.restLength, {
    this.compliance = 0.0,
  });

  final int nodeA;
  final int nodeB;
  final double restLength;
  final double compliance;
  double lambda = 0.0;

  void solve(List<VerletNode> nodes, double invDtSq) {
    lambda = _solveDistancePair(
      nodes,
      nodeA,
      nodeB,
      restLength,
      compliance,
      lambda,
      invDtSq,
    );
  }

  double currentError(List<VerletNode> nodes) =>
      (nodes[nodeA].position.distanceTo(nodes[nodeB].position) - restLength)
          .abs();
}

/// XPBD bend / structural cross-link constraint between two particle indices.
class VerletBendLink {
  VerletBendLink(
    this.nodeA,
    this.nodeB,
    this.restLength, {
    this.compliance = 0.02,
  });

  final int nodeA;
  final int nodeB;
  final double restLength;
  final double compliance;
  double lambda = 0.0;

  void solve(List<VerletNode> nodes, double invDtSq) {
    lambda = _solveDistancePair(
      nodes,
      nodeA,
      nodeB,
      restLength,
      compliance,
      lambda,
      invDtSq,
    );
  }

  double currentError(List<VerletNode> nodes) =>
      (nodes[nodeA].position.distanceTo(nodes[nodeB].position) - restLength)
          .abs();
}

/// XPBD signed-volume tetrahedral constraint across 4 particle indices.
class VerletTetVolume {
  VerletTetVolume(
    this.a,
    this.b,
    this.c,
    this.d,
    this.restVolume, {
    this.compliance = 0.001,
  });

  final int a;
  final int b;
  final int c;
  final int d;
  final double restVolume;
  final double compliance;
  double lambda = 0.0;

  /// Signed volume of tetrahedron (a, b, c, d):
  /// `(1 / 6) * (a - d) · ((b - d) × (c - d))`.
  static double signedTetVolume(
    VerletVec3 a,
    VerletVec3 b,
    VerletVec3 c,
    VerletVec3 d,
  ) {
    final ad = a - d;
    final bd = b - d;
    final cd = c - d;
    return ad.dot(bd.cross(cd)) / 6.0;
  }

  void solve(List<VerletNode> nodes, double invDtSq) {
    final nA = nodes[a];
    final nB = nodes[b];
    final nC = nodes[c];
    final nD = nodes[d];

    final wA = nA.isPinned ? 0.0 : nA.invMass;
    final wB = nB.isPinned ? 0.0 : nB.invMass;
    final wC = nC.isPinned ? 0.0 : nC.invMass;
    final wD = nD.isPinned ? 0.0 : nD.invMass;

    final pA = nA.position;
    final pB = nB.position;
    final pC = nC.position;
    final pD = nD.position;

    final ad = pA - pD;
    final bd = pB - pD;
    final cd = pC - pD;

    const inv6 = 1.0 / 6.0;
    final g0 = bd.cross(cd) * inv6;
    final g1 = cd.cross(ad) * inv6;
    final g2 = ad.cross(bd) * inv6;
    final g3 = -(g0 + g1 + g2);

    final wSum = wA * g0.lengthSquared +
        wB * g1.lengthSquared +
        wC * g2.lengthSquared +
        wD * g3.lengthSquared;

    final alphaTilde = compliance * invDtSq;
    final denom = wSum + alphaTilde;
    if (denom <= 1e-24) return;

    final volume = ad.dot(g0);
    final constraint = volume - restVolume;
    final dLambda = (-constraint - alphaTilde * lambda) / denom;
    lambda += dLambda;

    if (wA > 0.0) nA.position = pA + g0 * (wA * dLambda);
    if (wB > 0.0) nB.position = pB + g1 * (wB * dLambda);
    if (wC > 0.0) nC.position = pC + g2 * (wC * dLambda);
    if (wD > 0.0) nD.position = pD + g3 * (wD * dLambda);
  }

  double currentError(List<VerletNode> nodes) => (signedTetVolume(
            nodes[a].position,
            nodes[b].position,
            nodes[c].position,
            nodes[d].position,
          ) -
          restVolume)
      .abs();
}

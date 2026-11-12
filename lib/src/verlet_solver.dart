import 'dart:math' as math;

import 'verlet_constraints.dart';
import 'verlet_math.dart';

/// Substepped XPBD soft-body & Verlet constraint solver.
class VerletSoftBodySolver {
  VerletSoftBodySolver({
    List<VerletNode>? nodes,
    List<VerletDistanceLink>? distanceLinks,
    List<VerletBendLink>? bendLinks,
    List<VerletTetVolume>? tetVolumes,
    this.gravity = const VerletVec3(0, -9.81, 0),
    this.floorY = -1.0,
    this.damping = 0.995,
    this.substeps = 8,
    this.restitution = 0.2,
  })  : nodes = nodes ?? <VerletNode>[],
        distanceLinks = distanceLinks ?? <VerletDistanceLink>[],
        bendLinks = bendLinks ?? <VerletBendLink>[],
        tetVolumes = tetVolumes ?? <VerletTetVolume>[];

  final List<VerletNode> nodes;
  final List<VerletDistanceLink> distanceLinks;
  final List<VerletBendLink> bendLinks;
  final List<VerletTetVolume> tetVolumes;

  VerletVec3 gravity;
  double floorY;
  double damping;
  int substeps;
  double restitution;

  /// Builds an 8-vertex cube with 12 edges + 6 face-diagonals / structural
  /// links + 5 positive-volume tetrahedra covering the cube cleanly.
  static VerletSoftBodySolver createCubeRig({
    VerletVec3 center = VerletVec3.zero,
    double halfExtent = 0.5,
    double distanceCompliance = 0.001,
    double volumeCompliance = 0.0005,
    VerletVec3 gravity = const VerletVec3(0, -9.81, 0),
    double floorY = -1.0,
    int substeps = 8,
  }) {
    final h = halfExtent;
    final nodes = <VerletNode>[
      VerletNode(position: center + VerletVec3(-h, -h, -h)), // 0
      VerletNode(position: center + VerletVec3(h, -h, -h)), // 1
      VerletNode(position: center + VerletVec3(-h, h, -h)), // 2
      VerletNode(position: center + VerletVec3(h, h, -h)), // 3
      VerletNode(position: center + VerletVec3(-h, -h, h)), // 4
      VerletNode(position: center + VerletVec3(h, -h, h)), // 5
      VerletNode(position: center + VerletVec3(-h, h, h)), // 6
      VerletNode(position: center + VerletVec3(h, h, h)), // 7
    ];

    // 12 cube edges + 6 face diagonals belonging to the central tetrahedron.
    const distancePairs = <(int, int)>[
      // 12 edges
      (0, 1), (2, 3), (4, 5), (6, 7),
      (0, 2), (1, 3), (4, 6), (5, 7),
      (0, 4), (1, 5), (2, 6), (3, 7),
      // 6 face diagonals
      (0, 3), (5, 6), (0, 5), (3, 6), (0, 6), (3, 5),
    ];

    final distanceLinks = <VerletDistanceLink>[
      for (final (a, b) in distancePairs)
        VerletDistanceLink(
          a,
          b,
          nodes[a].position.distanceTo(nodes[b].position),
          compliance: distanceCompliance,
        ),
    ];

    // Opposite 6 face diagonals as symmetric shear/bend links.
    const bendPairs = <(int, int)>[
      (1, 2),
      (4, 7),
      (1, 4),
      (2, 7),
      (2, 4),
      (1, 7),
    ];

    final bendLinks = <VerletBendLink>[
      for (final (a, b) in bendPairs)
        VerletBendLink(
          a,
          b,
          nodes[a].position.distanceTo(nodes[b].position),
          compliance: distanceCompliance * 2.0,
        ),
    ];

    // 5-tetrahedron decomposition with positive signed volumes.
    const tetIndices = <(int, int, int, int)>[
      (0, 3, 5, 6),
      (1, 0, 3, 5),
      (2, 3, 0, 6),
      (4, 0, 5, 6),
      (7, 3, 6, 5),
    ];

    final tetVolumes = <VerletTetVolume>[
      for (final (a, b, c, d) in tetIndices)
        VerletTetVolume(
          a,
          b,
          c,
          d,
          VerletTetVolume.signedTetVolume(
            nodes[a].position,
            nodes[b].position,
            nodes[c].position,
            nodes[d].position,
          ),
          compliance: volumeCompliance,
        ),
    ];

    return VerletSoftBodySolver(
      nodes: nodes,
      distanceLinks: distanceLinks,
      bendLinks: bendLinks,
      tetVolumes: tetVolumes,
      gravity: gravity,
      floorY: floorY,
      substeps: substeps,
    );
  }

  /// Advances the XPBD simulation by [dt] seconds.
  void step(double dt) {
    if (dt <= 0.0 || !dt.isFinite) return;

    final clampedSubsteps = substeps < 1 ? 1 : substeps;
    final subDt = dt / clampedSubsteps;
    final invSubDt = 1.0 / subDt;
    final invDtSq = invSubDt * invSubDt;
    final clampedDamping = damping.clamp(0.0, 1.0);
    final subDamping = clampedDamping == 1.0
        ? 1.0
        : math.pow(clampedDamping, 1.0 / clampedSubsteps).toDouble();
    final clampedRestitution = restitution.clamp(0.0, 1.0);

    for (var s = 0; s < clampedSubsteps; s++) {
      // 1. Integrate velocity and predict positions.
      for (final node in nodes) {
        if (node.isPinned) {
          node.prevPosition = node.position;
          node.velocity = VerletVec3.zero;
          continue;
        }
        node.prevPosition = node.position;
        node.velocity = (node.velocity + gravity * subDt) * subDamping;
        node.position = node.position + node.velocity * subDt;
      }

      // 2. Clear Lagrange multipliers for this substep.
      for (final link in distanceLinks) {
        link.lambda = 0.0;
      }
      for (final link in bendLinks) {
        link.lambda = 0.0;
      }
      for (final tet in tetVolumes) {
        tet.lambda = 0.0;
      }

      // 3. Project XPBD constraints.
      for (final link in distanceLinks) {
        link.solve(nodes, invDtSq);
      }
      for (final link in bendLinks) {
        link.solve(nodes, invDtSq);
      }
      for (final tet in tetVolumes) {
        tet.solve(nodes, invDtSq);
      }

      // 4. Collide with floorY (with velocity restitution) and update velocity.
      for (final node in nodes) {
        if (node.isPinned) {
          node.velocity = VerletVec3.zero;
          continue;
        }
        var vel = (node.position - node.prevPosition) * invSubDt;
        if (node.position.y < floorY) {
          node.position = VerletVec3(node.position.x, floorY, node.position.z);
          if (vel.y < 0.0) {
            vel = VerletVec3(vel.x, -vel.y * clampedRestitution, vel.z);
          }
        }
        node.velocity = vel;
      }
    }
  }

  /// Applies an impulse to all unpinned nodes (`velocity += impulse * invMass`).
  void applyImpulse(VerletVec3 impulse) {
    for (final node in nodes) {
      if (!node.isPinned) {
        node.velocity = node.velocity + impulse * node.invMass;
      }
    }
  }

  /// Sum of absolute errors across all distance, bend, and volume constraints.
  double totalConstraintError() {
    var error = 0.0;
    for (final link in distanceLinks) {
      error += link.currentError(nodes);
    }
    for (final link in bendLinks) {
      error += link.currentError(nodes);
    }
    for (final tet in tetVolumes) {
      error += tet.currentError(nodes);
    }
    return error;
  }

  /// Mean normalized dimensionless strain across all active constraints.
  double averageStrain() {
    var strainSum = 0.0;
    var count = 0;

    for (final link in distanceLinks) {
      if (link.restLength > 1e-12) {
        strainSum += link.currentError(nodes) / link.restLength;
        count++;
      }
    }
    for (final link in bendLinks) {
      if (link.restLength > 1e-12) {
        strainSum += link.currentError(nodes) / link.restLength;
        count++;
      }
    }
    for (final tet in tetVolumes) {
      final denom = tet.restVolume.abs();
      if (denom > 1e-12) {
        strainSum += tet.currentError(nodes) / denom;
        count++;
      }
    }

    if (count == 0) return 0.0;
    return strainSum / count;
  }

  /// Mass-weighted center of mass of unpinned nodes (or geometric mean if all pinned).
  VerletVec3 centerOfMass() {
    if (nodes.isEmpty) return VerletVec3.zero;

    var weightedSum = VerletVec3.zero;
    var totalMass = 0.0;
    var geometricSum = VerletVec3.zero;

    for (final node in nodes) {
      geometricSum = geometricSum + node.position;
      if (!node.isPinned) {
        final mass = 1.0 / node.invMass;
        weightedSum = weightedSum + node.position * mass;
        totalMass += mass;
      }
    }

    if (totalMass <= 0.0) {
      return geometricSum * (1.0 / nodes.length);
    }
    return weightedSum * (1.0 / totalMass);
  }
}

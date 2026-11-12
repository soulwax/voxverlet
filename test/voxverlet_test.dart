import 'package:test/test.dart';
import 'package:voxverlet/voxverlet.dart';

void main() {
  group('VerletVec3 and VerletNode', () {
    test('vector algebra and invariants hold', () {
      const a = VerletVec3(1.0, 2.0, 3.0);
      const b = VerletVec3(4.0, -5.0, 6.0);

      expect(VerletVec3.zero.length, 0.0);
      expect(VerletVec3.zero.normalized, VerletVec3.zero);
      expect(a + b, const VerletVec3(5.0, -3.0, 9.0));
      expect(a - b, const VerletVec3(-3.0, 7.0, -3.0));
      expect(a * 2.0, const VerletVec3(2.0, 4.0, 6.0));
      expect(-a, const VerletVec3(-1.0, -2.0, -3.0));
      expect(a.dot(b), closeTo(12.0, 1e-12));
      expect(
        const VerletVec3(1, 0, 0).cross(const VerletVec3(0, 1, 0)),
        const VerletVec3(0, 0, 1),
      );
      expect(a.lengthSquared, closeTo(14.0, 1e-12));
      expect(a.normalized.length, closeTo(1.0, 1e-12));
      expect(a.distanceTo(a), 0.0);

      final pinned = VerletNode(position: a, invMass: 0.0);
      expect(pinned.isPinned, isTrue);
      final free = VerletNode(position: a);
      expect(free.isPinned, isFalse);
      expect(free.prevPosition, a);
      expect(free.velocity, VerletVec3.zero);
    });
  });

  group('VerletDistanceLink', () {
    test('converges stretched and compressed nodes toward restLength', () {
      const rest = 1.0;
      const invDtSq = 3600.0; // dt = 1/60

      // Stretched pair (initial distance = 2.0)
      final stretchedNodes = <VerletNode>[
        VerletNode(position: const VerletVec3(0, 0, 0)),
        VerletNode(position: const VerletVec3(2.0, 0, 0)),
      ];
      final stretchedLink = VerletDistanceLink(0, 1, rest, compliance: 0.0);
      stretchedLink.solve(stretchedNodes, invDtSq);

      expect(
        stretchedNodes[0].position.distanceTo(stretchedNodes[1].position),
        closeTo(rest, 1e-9),
      );
      // Center of mass preserved when equal masses.
      expect(
        (stretchedNodes[0].position.x + stretchedNodes[1].position.x) * 0.5,
        closeTo(1.0, 1e-9),
      );

      // Compressed pair (initial distance = 0.25)
      final compressedNodes = <VerletNode>[
        VerletNode(position: const VerletVec3(0, 0, 0)),
        VerletNode(position: const VerletVec3(0.25, 0, 0)),
      ];
      final compressedLink = VerletDistanceLink(0, 1, rest, compliance: 0.0);
      compressedLink.solve(compressedNodes, invDtSq);

      expect(
        compressedNodes[0].position.distanceTo(compressedNodes[1].position),
        closeTo(rest, 1e-9),
      );
    });

    test('respects compliance and pinned node constraints under load', () {
      // Top node pinned at origin, bottom node pulled down by gravity.
      VerletSoftBodySolver makeHangingLink(double compliance) {
        final nodes = <VerletNode>[
          VerletNode(position: VerletVec3.zero, invMass: 0.0),
          VerletNode(position: const VerletVec3(0, -1.0, 0), invMass: 1.0),
        ];
        final links = <VerletDistanceLink>[
          VerletDistanceLink(0, 1, 1.0, compliance: compliance),
        ];
        return VerletSoftBodySolver(
          nodes: nodes,
          distanceLinks: links,
          floorY: -100.0,
          substeps: 8,
        );
      }

      final stiff = makeHangingLink(0.0);
      final compliant = makeHangingLink(0.05);

      for (var i = 0; i < 60; i++) {
        stiff.step(1.0 / 60.0);
        compliant.step(1.0 / 60.0);
      }

      // Pinned anchor must not move.
      expect(stiff.nodes[0].position, VerletVec3.zero);
      expect(compliant.nodes[0].position, VerletVec3.zero);

      final stiffDist = stiff.nodes[0].position.distanceTo(
        stiff.nodes[1].position,
      );
      final compliantDist = compliant.nodes[0].position.distanceTo(
        compliant.nodes[1].position,
      );

      expect(stiffDist, closeTo(1.0, 1e-3));
      expect(compliantDist, greaterThan(stiffDist + 0.05));
    });
  });

  group('VerletTetVolume', () {
    test('preserves positive signed cube volume after compression and recovers on step', () {
      final rig = VerletSoftBodySolver.createCubeRig(
        center: const VerletVec3(0, 2.0, 0),
        halfExtent: 0.5,
        distanceCompliance: 0.0005,
        volumeCompliance: 0.0001,
        gravity: VerletVec3.zero,
        floorY: -10.0,
        substeps: 8,
      );

      expect(rig.nodes.length, 8);
      expect(rig.distanceLinks.length, 18);
      expect(rig.tetVolumes.length, 5);

      // Every tetrahedron must start with strictly positive signed volume,
      // and the 5 tetrahedra must sum to the full cube volume (1.0).
      var restTotalVolume = 0.0;
      for (final tet in rig.tetVolumes) {
        expect(tet.restVolume, greaterThan(0.0));
        restTotalVolume += tet.restVolume;
      }
      expect(restTotalVolume, closeTo(1.0, 1e-9));

      // Compress the cube vertically to 40% height around y = 2.0.
      for (final node in rig.nodes) {
        final dy = node.position.y - 2.0;
        node.position = VerletVec3(
          node.position.x,
          2.0 + dy * 0.4,
          node.position.z,
        );
        node.prevPosition = node.position;
      }

      double currentTotalVolume() {
        var sum = 0.0;
        for (final tet in rig.tetVolumes) {
          sum += VerletTetVolume.signedTetVolume(
            rig.nodes[tet.a].position,
            rig.nodes[tet.b].position,
            rig.nodes[tet.c].position,
            rig.nodes[tet.d].position,
          );
        }
        return sum;
      }

      final compressedVolume = currentTotalVolume();
      expect(compressedVolume, closeTo(0.4, 1e-9));

      // Step the solver; volume constraint should drive total volume back to 1.0
      // with every tetrahedron remaining positive.
      for (var i = 0; i < 40; i++) {
        rig.step(1.0 / 60.0);
      }

      final recoveredVolume = currentTotalVolume();
      expect(recoveredVolume, closeTo(1.0, 0.02));
      for (final tet in rig.tetVolumes) {
        final v = VerletTetVolume.signedTetVolume(
          rig.nodes[tet.a].position,
          rig.nodes[tet.b].position,
          rig.nodes[tet.c].position,
          rig.nodes[tet.d].position,
        );
        expect(v, greaterThan(0.0));
      }
      expect(rig.averageStrain(), lessThan(0.05));
    });
  });

  group('VerletSoftBodySolver.createCubeRig', () {
    test('falls under gravity, collides with floorY without NaNs/tunneling across 120 steps', () {
      const floorY = -1.0;
      final rig = VerletSoftBodySolver.createCubeRig(
        center: const VerletVec3(0, 1.5, 0),
        halfExtent: 0.5,
        floorY: floorY,
        substeps: 8,
      );

      final initialComY = rig.centerOfMass().y;
      expect(initialComY, closeTo(1.5, 1e-9));

      // Guard dt <= 0 does not alter state.
      rig.step(0.0);
      rig.step(-0.016);
      expect(rig.centerOfMass().y, closeTo(1.5, 1e-9));

      for (var i = 0; i < 120; i++) {
        rig.step(1.0 / 60.0);
        for (final node in rig.nodes) {
          expect(node.position.isFinite, isTrue);
          expect(node.velocity.isFinite, isTrue);
          expect(node.position.y, greaterThanOrEqualTo(floorY - 1e-9));
        }
      }

      final finalCom = rig.centerOfMass();
      expect(finalCom.isFinite, isTrue);
      // Cube should have fallen onto floorY = -1.0 (bottom at -1.0, COM around -0.5).
      expect(finalCom.y, lessThan(initialComY - 1.0));
      expect(finalCom.y, closeTo(floorY + 0.5, 0.15));
      expect(rig.totalConstraintError(), lessThan(0.25));
      expect(rig.averageStrain(), lessThan(0.05));

      // Impulse shifts velocity of unpinned nodes cleanly.
      rig.applyImpulse(const VerletVec3(2.0, 3.0, 0.0));
      for (final node in rig.nodes) {
        expect(node.velocity.x, greaterThan(1.5));
        expect(node.velocity.y, greaterThan(2.0));
      }
    });

    test('VerletBendLink reduces shear deformation', () {
      final nodes = <VerletNode>[
        VerletNode(position: const VerletVec3(0, 0, 0), invMass: 0.0),
        VerletNode(position: const VerletVec3(1, 0, 0)),
        VerletNode(position: const VerletVec3(1, 0.5, 0)),
      ];
      final bend = VerletBendLink(0, 2, 2.0, compliance: 0.0);
      bend.solve(nodes, 3600.0);
      expect(nodes[0].position.distanceTo(nodes[2].position), closeTo(2.0, 1e-9));
      expect(nodes[0].position, const VerletVec3(0, 0, 0));
    });
  });
}

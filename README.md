<p align="center">
  <img src="assets/logo.svg" width="88" alt="voxverlet logo"/>
</p>

# voxverlet — Pure-Dart XPBD soft-body & Verlet constraint solver for glVolition.

Self-contained, zero-dependency pure-Dart Extended Position Based Dynamics (XPBD)
soft-body and Verlet constraint solver package.

## Features

- `VerletVec3`: Immutable 3D vector math with zero allocations for standard ops.
- `VerletNode`: Particle state with position, previous position, velocity, inverse mass, and pinning.
- `VerletDistanceLink` & `VerletBendLink`: Compliant XPBD distance and structural/bend links.
- `VerletTetVolume`: Signed-volume tetrahedral XPBD constraint with exact analytical gradients.
- `VerletSoftBodySolver`: Substepped XPBD integrator with floor collision, restitution, impulse application, strain/error diagnostics, and `createCubeRig` factory.

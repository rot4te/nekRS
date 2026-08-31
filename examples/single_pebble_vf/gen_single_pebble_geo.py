#!/usr/bin/env python3
"""Generate a hex O-grid mesh (.geo) for flow around a single spherical
pebble centered in a square duct (paper Section 4.4).

Two-layer "cubed sphere" O-grid, both layers using the same 6-block
shared-edge topology (8 corner points, 12 edges, 6 face patches, 12 "fin"
side patches, 6 volumes) validated in concentric_spheres.geo:

  Layer 1 (curved): pebble surface (sphere, radius R) <-> a reference cube
  circumscribing it (half-side R, face centers touching the sphere).
  Corner edges are great-circle arcs (sphere side) connected by straight
  radial lines to the reference cube's straight edges.

  Layer 2 (flat): the same reference cube <-> the duct's outer walls (a
  rectangular box, half-extents Lx,Ly,Lz != R in general). Purely straight
  edges -- a plain nested-box extension, reusing Layer 1's outer corners as
  its inner corners.

Physical Surfaces: "pebble" (sphere), "inlet" (duct z=-Lz face), "outlet"
(duct z=+Lz face), "duct_wall" (the other 4 duct faces).
"""
import math

R  = 0.02   # pebble radius, m (diameter 4 cm)
RC = 0.02   # reference-cube half-side, m. Tried RC<R (e.g. 0.015) to avoid
            # the face-center "pinch" (zero shell thickness where the cube
            # face touches the sphere at RC=R) expecting fewer non-right-hand
            # elements from gmsh2nek's orientation check -- it produced MORE
            # (504 vs 96), so the pinch isn't the actual cause (a separate,
            # not-fully-root-caused Transfinite Volume orientation quirk is).
            # RC=R's 96 non-right-hand elements are fixed by gmsh2nek's
            # standard automatic correction and the resulting mesh validates
            # cleanly (0.025% view-factor closure error), so kept as is.
Lx = 0.04   # duct half-width, m (8 cm cross-section)
Ly = 0.04
Lz = 0.08   # duct half-length, m (16 cm length)

N_ARC  = 7   # divisions along each sphere/reference-cube edge (layer 1)
N_RAD1 = 5   # divisions through layer 1 (sphere -> reference cube)
N_RAD2 = 8   # divisions through layer 2 (reference cube -> duct wall)

signs = [
    (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
    (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
]
dirs = []
for (sx, sy, sz) in signs:
    n = math.sqrt(sx*sx + sy*sy + sz*sz)
    dirs.append((sx/n, sy/n, sz/n))

edges = [
    (1, 2), (2, 3), (3, 4), (4, 1),      # bottom
    (5, 6), (6, 7), (7, 8), (8, 5),      # top
    (1, 5), (2, 6), (3, 7), (4, 8),      # vertical
]
faces = {
    "bottom": [1, 2, 3, 4],
    "top":    [5, 6, 7, 8],
    "front":  [1, 10, -5, -9],
    "back":   [-3, 11, 7, -12],
    "left":   [-4, 12, 8, -9],
    "right":  [2, 11, -6, -10],
}
face_names = list(faces.keys())

lines = []
def emit(s): lines.append(s)

emit("// Single-pebble Nek-VF case (paper Section 4.4): square duct with a")
emit("// spherical pebble at the center. Two-layer cubed-sphere/cubed-box")
emit("// O-grid -- see the header of the generating script for the topology.")
emit("Mesh.SubdivisionAlgorithm = 0;")
emit(f"N_ARC = {N_ARC};  N_RAD1 = {N_RAD1};  N_RAD2 = {N_RAD2};")
emit("")
emit("Point(1000) = {0, 0, 0};  // origin")

# --- Layer 1 inner (sphere) points 101..108, outer (ref cube) points 201..208
for i, (dx, dy, dz) in enumerate(dirs, start=1):
    emit(f"Point({100+i}) = {{{dx*R:.10f}, {dy*R:.10f}, {dz*R:.10f}}};")
for i, (sx, sy, sz) in enumerate(signs, start=1):
    emit(f"Point({200+i}) = {{{sx*RC:.10f}, {sy*RC:.10f}, {sz*RC:.10f}}};")
# --- Layer 2 outer (duct box) points 301..308 -- same directions, scaled
#     per-axis to the duct's actual (Lx,Ly,Lz) half-extents.
for i, (sx, sy, sz) in enumerate(signs, start=1):
    emit(f"Point({300+i}) = {{{sx*Lx:.10f}, {sy*Ly:.10f}, {sz*Lz:.10f}}};")

emit("")
emit("// Layer 1 radial spokes: sphere corner -> reference-cube corner")
for i in range(1, 9):
    emit(f"Line({400+i}) = {{{100+i}, {200+i}}};")
emit("// Layer 2 radial spokes: reference-cube corner -> duct-box corner")
for i in range(1, 9):
    emit(f"Line({500+i}) = {{{200+i}, {300+i}}};")

emit("")
emit("// Sphere edges (great-circle arcs) and reference-cube edges (straight)")
for idx, (a, b) in enumerate(edges, start=1):
    emit(f"Circle({600+idx}) = {{{100+a}, 1000, {100+b}}};")
for idx, (a, b) in enumerate(edges, start=1):
    emit(f"Line({700+idx}) = {{{200+a}, {200+b}}};")
emit("// Duct-box edges (straight)")
for idx, (a, b) in enumerate(edges, start=1):
    emit(f"Line({800+idx}) = {{{300+a}, {300+b}}};")

emit("")
emit(f"Transfinite Curve{{{', '.join(str(600+i) for i in range(1,13))}, "
     f"{', '.join(str(700+i) for i in range(1,13))}}} = N_ARC;")
emit(f"Transfinite Curve{{{', '.join(str(800+i) for i in range(1,13))}}} = N_ARC;")
emit(f"Transfinite Curve{{{', '.join(str(400+i) for i in range(1,9))}}} = N_RAD1;")
emit(f"Transfinite Curve{{{', '.join(str(500+i) for i in range(1,9))}}} = N_RAD2;")

def curve_id(base, e):
    return (base + abs(e)) if e > 0 else -(base + abs(e))

emit("")
emit("// Layer 1 face patches: sphere (curved) and reference cube (flat)")
sphere_surf, refcube_surf = {}, {}
for fi, name in enumerate(face_names, start=1):
    eids = faces[name]
    loop_sph = ", ".join(str(curve_id(600, e)) for e in eids)
    loop_ref = ", ".join(str(curve_id(700, e)) for e in eids)
    lid_s, lid_r = 900+fi, 1000+fi
    sid_s, sid_r = 1100+fi, 1200+fi
    emit(f"Line Loop({lid_s}) = {{{loop_sph}}};")
    emit(f"Ruled Surface({sid_s}) = {{{lid_s}}};")
    emit(f"Line Loop({lid_r}) = {{{loop_ref}}};")
    emit(f"Ruled Surface({sid_r}) = {{{lid_r}}};")
    for sid in (sid_s, sid_r):
        emit(f"Transfinite Surface{{{sid}}};  Recombine Surface{{{sid}}};")
    sphere_surf[name], refcube_surf[name] = sid_s, sid_r

emit("")
emit("// Layer 1 fin (side) surfaces, one per cube edge, shared between the")
emit("// two face-blocks meeting at that edge.")
fin1 = {}
for idx, (a, b) in enumerate(edges, start=1):
    lid, sid = 1300+idx, 1400+idx
    emit(f"Line Loop({lid}) = {{{600+idx}, {400+b}, -{700+idx}, -{400+a}}};")
    emit(f"Ruled Surface({sid}) = {{{lid}}};")
    emit(f"Transfinite Surface{{{sid}}};  Recombine Surface{{{sid}}};")
    fin1[idx] = sid

emit("")
emit("// Layer 1 volumes (sphere shell)")
vol1 = {}
for fi, name in enumerate(face_names, start=1):
    eids = faces[name]
    sides = [fin1[abs(e)] for e in eids]
    surfs = [sphere_surf[name], refcube_surf[name]] + sides
    lid, vid = 1500+fi, 1600+fi
    emit(f"Surface Loop({lid}) = {{{', '.join(str(s) for s in surfs)}}};")
    emit(f"Volume({vid}) = {{{lid}}};")
    emit(f"Transfinite Volume{{{vid}}};  Recombine Volume{{{vid}}};")
    vol1[name] = vid

emit("")
emit("// Layer 2 duct-box face patches (flat)")
ductbox_surf = {}
for fi, name in enumerate(face_names, start=1):
    eids = faces[name]
    loop = ", ".join(str(curve_id(800, e)) for e in eids)
    lid, sid = 1700+fi, 1800+fi
    emit(f"Line Loop({lid}) = {{{loop}}};")
    emit(f"Ruled Surface({sid}) = {{{lid}}};")
    emit(f"Transfinite Surface{{{sid}}};  Recombine Surface{{{sid}}};")
    ductbox_surf[name] = sid

emit("")
emit("// Layer 2 fin surfaces (reference-cube edge -> duct-box edge)")
fin2 = {}
for idx, (a, b) in enumerate(edges, start=1):
    lid, sid = 1900+idx, 2000+idx
    emit(f"Line Loop({lid}) = {{{700+idx}, {500+b}, -{800+idx}, -{500+a}}};")
    emit(f"Ruled Surface({sid}) = {{{lid}}};")
    emit(f"Transfinite Surface{{{sid}}};  Recombine Surface{{{sid}}};")
    fin2[idx] = sid

emit("")
emit("// Layer 2 volumes (reference cube -> duct wall)")
vol2 = {}
for fi, name in enumerate(face_names, start=1):
    eids = faces[name]
    sides = [fin2[abs(e)] for e in eids]
    surfs = [refcube_surf[name], ductbox_surf[name]] + sides
    lid, vid = 2100+fi, 2200+fi
    emit(f"Surface Loop({lid}) = {{{', '.join(str(s) for s in surfs)}}};")
    emit(f"Volume({vid}) = {{{lid}}};")
    emit(f"Transfinite Volume{{{vid}}};  Recombine Volume{{{vid}}};")
    vol2[name] = vid

emit("")
emit(f"Physical Surface(\"pebble\", 1) = {{{', '.join(str(sphere_surf[n]) for n in face_names)}}};")
emit(f"Physical Surface(\"inlet\", 2) = {{{ductbox_surf['bottom']}}};")   # z=-Lz
emit(f"Physical Surface(\"outlet\", 3) = {{{ductbox_surf['top']}}};")     # z=+Lz
emit("Physical Surface(\"duct_wall\", 4) = {" +
     ", ".join(str(ductbox_surf[n]) for n in ("front","back","left","right")) + "};")
emit("Physical Volume(\"fluid\", 5) = {" +
     ", ".join(str(vol1[n]) for n in face_names) + ", " +
     ", ".join(str(vol2[n]) for n in face_names) + "};")

with open("/Users/coxea3/nekRS/single_pebble_vf/single_pebble.geo", "w") as f:
    f.write("\n".join(lines) + "\n")
print("wrote single_pebble.geo")

#!/usr/bin/env python3
"""Generate a cubed-sphere hex mesh (.geo) for two concentric spherical
shells: inner radius r1 ("hot"), outer radius r2 ("cold"), with the
enclosed volume between them as the fluid/void domain. Six curved hex
blocks tile the shell, each bounded by an inner sphere patch, an outer
sphere patch, and four "radial fin" side patches shared with neighboring
blocks (standard cubed-sphere O-grid decomposition).
"""
import math

r1 = 0.5   # inner (hot) sphere radius, m -- derived from Table II's A1/A2=0.25
r2 = 1.0   # outer (cold) sphere radius, m
N_ARC = 16  # divisions along each great-circle arc edge
N_RAD = 4  # divisions in the radial direction

# 8 cube-corner directions (unit vectors), standard hex corner order:
# bottom face 1,2,3,4 (z<0) then top face 5,6,7,8 (z>0), 5 above 1 etc.
signs = [
    (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
    (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
]
dirs = []
for (sx, sy, sz) in signs:
    v = (sx, sy, sz)
    n = math.sqrt(sx*sx + sy*sy + sz*sz)
    dirs.append((v[0]/n, v[1]/n, v[2]/n))

# 12 cube edges as corner-index pairs (1-indexed corners, matching `signs`)
edges = [
    (1, 2), (2, 3), (3, 4), (4, 1),      # bottom
    (5, 6), (6, 7), (7, 8), (8, 5),      # top
    (1, 5), (2, 6), (3, 7), (4, 8),      # vertical
]

# 6 cube faces as ordered 4-corner loops (consistent winding), and the
# corresponding 4 edge indices (1-based into `edges`, signed by direction
# matching the corner-loop traversal)
faces = {
    "bottom": ([1, 2, 3, 4], [1, 2, 3, 4]),        # z=-r
    "top":    ([5, 6, 7, 8], [5, 6, 7, 8]),        # z=+r
    "front":  ([1, 2, 6, 5], [1, 10, -5, -9]),     # y=-r  (1-2,2-6,6-5,5-1)
    "back":   ([4, 3, 7, 8], [-3, 11, 7, -12]),    # y=+r  (4-3,3-7,7-8,8-4)
    "left":   ([1, 4, 8, 5], [-4, 12, 8, -9]),     # x=-r  (1-4,4-8,8-5,5-1)
    "right":  ([2, 3, 7, 6], [2, 11, -6, -10]),    # x=+r  (2-3,3-7,7-6,6-2)
}

lines = []
def emit(s):
    lines.append(s)

emit(f"// Concentric-spheres Nek-VF verification case (paper Section 4.2).")
emit(f"// Cubed-sphere shell mesh: 6 curved hex blocks between r1={r1} (hot)")
emit(f"// and r2={r2} (cold). r1/r2=0.5 is derived from Table II's implied")
emit(f"// A1/A2=0.25 (paper states neither radius explicitly).")
emit("Mesh.SubdivisionAlgorithm = 0;")
emit(f"N_ARC = {N_ARC};")
emit(f"N_RAD = {N_RAD};")
emit("")
emit("Point(100) = {0, 0, 0};  // sphere center")

# 8 inner points (101-108), 8 outer points (201-208)
for i, (dx, dy, dz) in enumerate(dirs, start=1):
    emit(f"Point({100+i}) = {{{dx*r1:.10f}, {dy*r1:.10f}, {dz*r1:.10f}}};")
for i, (dx, dy, dz) in enumerate(dirs, start=1):
    emit(f"Point({200+i}) = {{{dx*r2:.10f}, {dy*r2:.10f}, {dz*r2:.10f}}};")

emit("")
# 8 radial lines: 301..308, connecting inner corner i -> outer corner i
for i in range(1, 9):
    emit(f"Line({300+i}) = {{{100+i}, {200+i}}};")

emit("")
# 12 inner arcs (401..412) and 12 outer arcs (501..512), one per cube edge
for idx, (a, b) in enumerate(edges, start=1):
    emit(f"Circle({400+idx}) = {{{100+a}, 100, {100+b}}};")
for idx, (a, b) in enumerate(edges, start=1):
    emit(f"Circle({500+idx}) = {{{200+a}, 100, {200+b}}};")

emit("")
emit("Transfinite Curve{" + ", ".join(str(400+i) for i in range(1, 13)) +
     ", " + ", ".join(str(500+i) for i in range(1, 13)) + "} = N_ARC;")
emit("Transfinite Curve{" + ", ".join(str(300+i) for i in range(1, 9)) +
     "} = N_RAD;")

emit("")
emit("// Inner (hot) and outer (cold) sphere patches, one per cube face")
face_names = list(faces.keys())
inner_surf = {}
outer_surf = {}
for fi, name in enumerate(face_names, start=1):
    corners, edge_ids = faces[name]
    loop_in = ", ".join(str(400 + abs(e)) if e > 0 else str(-(400 + abs(e))) for e in edge_ids)
    loop_out = ", ".join(str(500 + abs(e)) if e > 0 else str(-(500 + abs(e))) for e in edge_ids)
    lid_in = 600 + fi
    lid_out = 700 + fi
    sid_in = 800 + fi
    sid_out = 900 + fi
    emit(f"Line Loop({lid_in}) = {{{loop_in}}};")
    emit(f"Ruled Surface({sid_in}) = {{{lid_in}}};")
    emit(f"Line Loop({lid_out}) = {{{loop_out}}};")
    emit(f"Ruled Surface({sid_out}) = {{{lid_out}}};")
    emit(f"Transfinite Surface{{{sid_in}}};")
    emit(f"Transfinite Surface{{{sid_out}}};")
    emit(f"Recombine Surface{{{sid_in}}};")
    emit(f"Recombine Surface{{{sid_out}}};")
    inner_surf[name] = sid_in
    outer_surf[name] = sid_out

emit("")
emit("// 12 radial 'fin' side surfaces, one per cube edge, shared between")
emit("// the two hex blocks (cube faces) that meet at that edge.")
side_surf = {}
for idx, (a, b) in enumerate(edges, start=1):
    # side surface bounded by: inner arc (a->b), radial line at b, outer
    # arc reversed (b->a), radial line at a reversed
    lid = 1000 + idx
    sid = 1100 + idx
    emit(f"Line Loop({lid}) = {{{400+idx}, {300+b}, -{500+idx}, -{300+a}}};")
    emit(f"Ruled Surface({sid}) = {{{lid}}};")
    emit(f"Transfinite Surface{{{sid}}};")
    emit(f"Recombine Surface{{{sid}}};")
    side_surf[idx] = sid

emit("")
emit("// 6 hex volumes, one per cube face, bounded by its inner patch,")
emit("// outer patch, and its 4 (shared) side patches.")
vol_ids = {}
for fi, name in enumerate(face_names, start=1):
    corners, edge_ids = faces[name]
    sides = [side_surf[abs(e)] for e in edge_ids]
    surfs = [inner_surf[name], outer_surf[name]] + sides
    lid = 1200 + fi
    vid = 1300 + fi
    emit(f"Surface Loop({lid}) = {{{', '.join(str(s) for s in surfs)}}};")
    emit(f"Volume({vid}) = {{{lid}}};")
    emit(f"Transfinite Volume{{{vid}}};")
    emit(f"Recombine Volume{{{vid}}};")
    vol_ids[name] = vid

emit("")
emit(f"Physical Surface(\"hot\", 1) = {{{', '.join(str(inner_surf[n]) for n in face_names)}}};")
emit(f"Physical Surface(\"cold\", 2) = {{{', '.join(str(outer_surf[n]) for n in face_names)}}};")
emit(f"Physical Volume(\"fluid\", 3) = {{{', '.join(str(vol_ids[n]) for n in face_names)}}};")

with open("/Users/coxea3/nekRS/concentric_spheres_vf/concentric_spheres.geo", "w") as f:
    f.write("\n".join(lines) + "\n")

print("wrote concentric_spheres.geo")
print("r1=", r1, "r2=", r2, "A1/A2=", (r1/r2)**2)

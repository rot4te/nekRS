#!/usr/bin/env python3
"""
Generates singlePebble.geo: a structured hex mesh of flow around a single
spherical pebble centered inside a square duct, for the Yuan et al. Sec 4.4
("Flow around a single pebble") radiation-convection benchmark.

genbox can't do curved geometry (the pebble is a sphere), so this follows the
same gmsh "cubed sphere" technique already validated in
examples/concentricSpheres/make_mesh.py, extended to TWO nested shells
instead of one:

  Shell A (inner): pebble surface (sphere, radius Rp) -> a small cubed
  "collar" (flat cube, half-width Rc). Same technique as concentricSpheres:
  6 curved panels, each bounded by 2 great-circle caps (sphere side) and a
  flat cap (collar side), with straight radial connectors between them.

  Shell B (outer): the collar cube (half-width Rc) -> the duct's own outer
  walls (a box, half-widths Wx=Wy (duct cross-section) and Wz (duct half
  length) -- NOT a cube, since the duct is longer than it is wide). This
  shell is fully flat/linear (both its inner and outer boundary are planar
  quads), so it's a simpler case of the same 6-panel construction: a nested
  box-in-box nested-shell mesh, sharing face/edge topology with Shell A's
  collar so the two shells are exactly conformal at the collar interface.

Physical Surface groups (raw boundary IDs, see singlePebble.par):
  1 = pebble surface (all 6 Shell-A sphere caps)
  2 = duct side walls (Shell-B duct caps for the +-x, +-y panels)
  3 = inlet (Shell-B duct cap for the -z panel)
  4 = outlet (Shell-B duct cap for the +z panel)
The 6 collar-cube caps and all 24 "side" surfaces (12 per shell) are left
untagged -- they are internal mesh interfaces (Shell A/Shell B share the
collar caps verbatim, by construction, not by a coincidence-merge step), not
domain boundaries.
"""
import itertools
import math

# ---- geometry parameters (meters, matching the paper's cm dimensions) ----
RP = 0.02          # pebble radius (D=4 cm, Table IV/Sec 4.4)
RC = 0.028          # collar cube half-width (must be > RP; chosen for a
                    # modest, not-too-distorted Shell A)
WX = 0.04           # duct half cross-section (8x8 cm^2 duct, Sec 4.4)
WY = 0.04
WZ = 0.08           # duct half-length (16 cm long, pebble centered at z=0)

# ---- mesh resolution ----
# N: points per panel edge (>=2), applied uniformly to EVERY tangential
# edge in both shells (sphere, collar, and duct caps all share edges
# pairwise across panels/shells, so a single global N keeps the whole
# structure conformal -- same constraint as concentricSpheres' single N).
# N=5 -> 4 elements/edge, 16 elements/panel cap, matching or exceeding
# concentricSpheres' N=4 curvature-error validation (0.22% at N=4).
N = 5
# NR_A: radial points, pebble surface -> collar (>=2). NR_A=3 -> 2 elements.
NR_A = 3
# NR_B: radial points, collar -> duct outer boundary (>=2). Uniform across
# all 6 Shell-B panels because each panel's radial edges are SHARED with
# its two neighboring panels at each corner (same O-grid closure constraint
# as concentricSpheres) -- independently grading the short x/y extension
# (1.2 cm) vs. the long z extension (5.2 cm) would break that sharing.
# NR_B=3 -> 2 elements; accepted as a coarse-mesh compromise for this demo
# case (documented in singlePebble.md), not a mesh-converged production mesh.
NR_B = 3

lc = 1.0  # unused (transfinite everywhere), required by Point syntax

lines = []
def emit(s):
    lines.append(s)

emit(f"// singlePebble: flow around a single pebble (Yuan et al. Sec 4.4)")
emit(f"// Rp={RP}, Rc={RC}, duct half-widths=({WX},{WY},{WZ})")
emit("SetFactory(\"Built-in\");")
emit("")

emit(f"Point(1) = {{0, 0, 0, {lc}}};")
ORIGIN = 1
next_pt = 2

signs = list(itertools.product([-1, 1], repeat=3))  # (sx,sy,sz)

# ---- 3 layers: sphere (curved), collar (flat cube), duct (flat box) ----
layers = [
    {"name": "sphere", "kind": "sphere", "R": RP},
    {"name": "collar", "kind": "box", "W": (RC, RC, RC)},
    {"name": "duct", "kind": "box", "W": (WX, WY, WZ)},
]

def corner_xyz(layer, s):
    sx, sy, sz = s
    if layer["kind"] == "sphere":
        R = layer["R"]
        sqrt3 = math.sqrt(3)
        return (sx * R / sqrt3, sy * R / sqrt3, sz * R / sqrt3)
    else:
        Wx, Wy, Wz = layer["W"]
        return (sx * Wx, sy * Wy, sz * Wz)

# corner points per layer
corner_pt = {}  # (layerIdx, s) -> point id
for li, layer in enumerate(layers):
    for s in signs:
        pid = next_pt
        next_pt += 1
        x, y, z = corner_xyz(layer, s)
        corner_pt[(li, s)] = pid
        emit(f"Point({pid}) = {{{x:.10f}, {y:.10f}, {z:.10f}, {lc}}};")
emit("")

# ---- 12 tangential edges per layer (arc if sphere, line if box) ----
def other_axes(axis):
    return [a for a in range(3) if a != axis]

def corner_from_fixed(varying_axis, varying_sign, fixed):
    s = [0, 0, 0]
    s[varying_axis] = varying_sign
    for a, sgn in fixed.items():
        s[a] = sgn
    return tuple(s)

edges = []  # (varying_axis, fixed_axis0, fixed_sign0, fixed_axis1, fixed_sign1)
for varying_axis in range(3):
    fa0, fa1 = other_axes(varying_axis)
    for s0, s1 in itertools.product([-1, 1], repeat=2):
        edges.append((varying_axis, fa0, s0, fa1, s1))

next_line = 1
layer_edge_arc = {}  # (layerIdx, edgeKey) -> (cA, cB, lineID)
for li, layer in enumerate(layers):
    for key in edges:
        va, fa0, s0, fa1, s1 = key
        fixed = {fa0: s0, fa1: s1}
        cA = corner_from_fixed(va, -1, fixed)
        cB = corner_from_fixed(va, +1, fixed)
        pA = corner_pt[(li, cA)]
        pB = corner_pt[(li, cB)]
        lid = next_line
        next_line += 1
        if layer["kind"] == "sphere":
            emit(f"Circle({lid}) = {{{pA}, {ORIGIN}, {pB}}};")
        else:
            emit(f"Line({lid}) = {{{pA}, {pB}}};")
        layer_edge_arc[(li, key)] = (cA, cB, lid)
        emit(f"Transfinite Curve {{{lid}}} = {N};")
emit("")

# ---- radial connector lines: layer li -> li+1, one per octant corner ----
# NR_for_shell[shellIdx] gives the point count for that shell's radial lines
NR_for_shell = [NR_A, NR_B]
shell_radial_line = [{}, {}]  # shellIdx -> {s: lineID}
for shell in (0, 1):
    li_in, li_out = shell, shell + 1
    for s in signs:
        pin = corner_pt[(li_in, s)]
        pout = corner_pt[(li_out, s)]
        lid = next_line
        next_line += 1
        emit(f"Line({lid}) = {{{pin}, {pout}}};")
        emit(f"Transfinite Curve {{{lid}}} = {NR_for_shell[shell]};")
        shell_radial_line[shell][s] = lid
emit("")

# ---- panel corner ordering (CCW as seen from outside), shared by all layers ----
panel_defs = []
for axis in range(3):
    for sign in (-1, 1):
        panel_defs.append((axis, sign))

def panel_corners_ccw(axis, sign):
    free = other_axes(axis)
    a0, a1 = free
    order = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    corners = []
    for (s0, s1) in order:
        s = [0, 0, 0]
        s[axis] = sign
        s[a0] = s0
        s[a1] = s1
        corners.append(tuple(s))
    return corners, a0, a1

def find_edge_key(cA, cB):
    diffs = [i for i in range(3) if cA[i] != cB[i]]
    assert len(diffs) == 1
    va = diffs[0]
    fa0, fa1 = other_axes(va)
    fixed = {fa0: cA[fa0], fa1: cA[fa1]}
    return (va, fa0, fixed[fa0], fa1, fixed[fa1])

# ---- cap surfaces: one per (layer, panel), built ONCE, shared between
# whichever shell(s) touch that layer (layer 0: shell0 inner only; layer 1:
# shell0 outer AND shell1 inner -- same surface; layer 2: shell1 outer only)
surf_lines = []
next_surf = 1
layer_panel_cap = {}  # (li, panelIdx) -> surfID

for li, layer in enumerate(layers):
    for pidx, (axis, sign) in enumerate(panel_defs):
        corners, a0, a1 = panel_corners_ccw(axis, sign)
        loop_lines = []
        for i in range(4):
            cA = corners[i]
            cB = corners[(i + 1) % 4]
            key = find_edge_key(cA, cB)
            cAk, cBk, lid = layer_edge_arc[(li, key)]
            loop_lines.append(lid if cA == cAk else -lid)
        sid = next_surf
        next_surf += 1
        cl_id = 10000 + sid
        surf_lines.append(f"Curve Loop({cl_id}) = {{{','.join(str(x) for x in loop_lines)}}};")
        surf_lines.append(f"Surface({sid}) = {{{cl_id}}};")
        pts_in_order = [corner_pt[(li, c)] for c in corners]
        surf_lines.append(f"Transfinite Surface {{{sid}}} = {{{','.join(str(p) for p in pts_in_order)}}};")
        surf_lines.append(f"Recombine Surface {{{sid}}};")
        layer_panel_cap[(li, pidx)] = sid

# ---- side surfaces: one per (shell, panel, of its 4 tangential edges),
# shared between adjacent panels WITHIN a shell (same corner-edge, differing
# panel), never shared across shells (radial band is shell-specific)
shell_side_surf = [{}, {}]  # shellIdx -> {edgeKey: surfID}
for shell in (0, 1):
    li_in, li_out = shell, shell + 1
    for pidx, (axis, sign) in enumerate(panel_defs):
        corners, a0, a1 = panel_corners_ccw(axis, sign)
        for i in range(4):
            cA = corners[i]
            cB = corners[(i + 1) % 4]
            key = find_edge_key(cA, cB)
            if key in shell_side_surf[shell]:
                continue
            _, _, lid_in = layer_edge_arc[(li_in, key)]
            cAk, cBk, lid_out = layer_edge_arc[(li_out, key)]
            lineA = shell_radial_line[shell][cAk]  # cAk_in -> cAk_out
            lineB = shell_radial_line[shell][cBk]
            sid = next_surf
            next_surf += 1
            cl_id = 10000 + sid
            surf_lines.append(f"Curve Loop({cl_id}) = {{{lid_in},{lineB},{-lid_out},{-lineA}}};")
            surf_lines.append(f"Surface({sid}) = {{{cl_id}}};")
            pin_A = corner_pt[(li_in, cAk)]
            pin_B = corner_pt[(li_in, cBk)]
            pout_B = corner_pt[(li_out, cBk)]
            pout_A = corner_pt[(li_out, cAk)]
            surf_lines.append(f"Transfinite Surface {{{sid}}} = {{{pin_A},{pin_B},{pout_B},{pout_A}}};")
            surf_lines.append(f"Recombine Surface {{{sid}}};")
            shell_side_surf[shell][key] = sid

emit("\n".join(surf_lines))
emit("")

# ---- volumes: one per (shell, panel) ----
vol_lines = []
shell_panel_vol = [{}, {}]
for shell in (0, 1):
    li_in, li_out = shell, shell + 1
    for pidx, (axis, sign) in enumerate(panel_defs):
        corners, a0, a1 = panel_corners_ccw(axis, sign)
        sid_in = layer_panel_cap[(li_in, pidx)]
        sid_out = layer_panel_cap[(li_out, pidx)]
        side_sids = []
        for i in range(4):
            cA = corners[i]
            cB = corners[(i + 1) % 4]
            key = find_edge_key(cA, cB)
            side_sids.append(shell_side_surf[shell][key])
        vid = 1000 * (shell + 1) + sid_in
        sl_id = 20000 + vid
        vol_lines.append(f"Surface Loop({sl_id}) = {{{sid_in},{sid_out},{','.join(str(s) for s in side_sids)}}};")
        vol_lines.append(f"Volume({vid}) = {{{sl_id}}};")
        pts_in = [corner_pt[(li_in, c)] for c in corners]
        pts_out = [corner_pt[(li_out, c)] for c in corners]
        vol_lines.append(f"Transfinite Volume {{{vid}}} = {{{','.join(str(p) for p in pts_in + pts_out)}}};")
        vol_lines.append(f"Recombine Volume {{{vid}}};")
        shell_panel_vol[shell][pidx] = vid

emit("\n".join(vol_lines))
emit("")

# ---- physical groups ----
pebble_surfs = [layer_panel_cap[(0, pidx)] for pidx in range(6)]

wall_pidx = []
inlet_pidx = None
outlet_pidx = None
for pidx, (axis, sign) in enumerate(panel_defs):
    if axis == 2 and sign == -1:
        inlet_pidx = pidx
    elif axis == 2 and sign == 1:
        outlet_pidx = pidx
    else:
        wall_pidx.append(pidx)

wall_surfs = [layer_panel_cap[(2, pidx)] for pidx in wall_pidx]
inlet_surf = layer_panel_cap[(2, inlet_pidx)]
outlet_surf = layer_panel_cap[(2, outlet_pidx)]

all_vols = list(shell_panel_vol[0].values()) + list(shell_panel_vol[1].values())

emit(f"Physical Surface(\"pebble\", 1) = {{{','.join(str(s) for s in pebble_surfs)}}};")
emit(f"Physical Surface(\"wall\", 2) = {{{','.join(str(s) for s in wall_surfs)}}};")
emit(f"Physical Surface(\"inlet\", 3) = {{{inlet_surf}}};")
emit(f"Physical Surface(\"outlet\", 4) = {{{outlet_surf}}};")
emit(f"Physical Volume(\"fluid\", 1) = {{{','.join(str(v) for v in all_vols)}}};")
emit("")

emit("Mesh.ElementOrder = 2;")
emit("Mesh.SecondOrderIncomplete = 1;")
emit("Mesh.RecombineAll = 1;")
# Override any ambient ~/.gmsh-options state (see concentricSpheres/make_mesh.py
# Problem #1) that would silently subdivide every transfinite element.
emit("Mesh.SubdivisionAlgorithm = 0;")
emit("Coherence;")

with open("singlePebble.geo", "w") as f:
    f.write("\n".join(lines) + "\n")

nElemA = 6 * (N - 1) * (N - 1) * (NR_A - 1)
nElemB = 6 * (N - 1) * (N - 1) * (NR_B - 1)
print("wrote singlePebble.geo")
print(f"points: {next_pt-1}, curves: {next_line-1}, surfaces: {next_surf-1}, volumes: {len(all_vols)}")
print(f"elements: Shell A (pebble->collar) = {nElemA}, Shell B (collar->duct) = {nElemB}, total = {nElemA+nElemB}")
print(f"pebble (bID=1) faces: {6*(N-1)*(N-1)}, wall (bID=2) faces: {4*(N-1)*(N-1)}, "
      f"inlet (bID=3) faces: {(N-1)*(N-1)}, outlet (bID=4) faces: {(N-1)*(N-1)}")

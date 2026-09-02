// Override a leftover global Gmsh GUI preference (~/.gmsh-options has
// Mesh.SubdivisionAlgorithm=2 from unrelated work) that would otherwise
// silently 8x-subdivide every transfinite hex on CLI invocation too.
Mesh.SubdivisionAlgorithm = 0;

// Parallel-plates Nek-VF verification case (paper Section 4.1).
// Closed box: hot plate (z=0), cold plate (z=H), four side walls.
// The paper gives wall temperatures and emissivity but not absolute plate
// dimensions or spacing -- Eq. 11 (infinite parallel plates) is only exact
// in the limit of plate size >> spacing, so a box with a large cross-section
// to height ratio is used here to approximate that limit while still being a
// genuine closed 6-surface enclosure that RadiativeViewFactor.jl computes
// exact (non-idealized) view factors for.
//
// W/H was originally 4 (W=4.0), which is NOT a good approximation of the
// idealized limit: side-wall total area (4*W*H) came out exactly equal to
// one plate's area (W*W) at that ratio, and a run against that mesh
// deviated from the paper's Table I by ~25% on the cold/side walls (vs
// ~5.5% on the hot wall). Raised to W/H=20 (2026-09-03) to shrink the
// side walls' relative contribution to 1/5 of a plate's area instead of
// 1x -- see nekRS/changelog.md (2026-09-03) for the run that motivated
// this. N_W/N_H left unchanged: this only changes the physical size of
// each element, not the element or DoF count (element count is set by the
// transfinite curve node counts below, independent of W/H) -- polynomial
// order (case.par) provides the actual within-element resolution, so
// there's no accuracy reason to also refine the mesh here.
W = 20.0;  // plate side length (m) -- assumed, not stated in the paper
H = 1.0;   // plate spacing (m)     -- assumed, not stated in the paper
N_W = 21;  // nodes across each plate edge
N_H = 11;  // nodes through the spacing direction

Point(1) = {0,0,0};
Point(2) = {W,0,0};
Point(3) = {W,W,0};
Point(4) = {0,W,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};
Transfinite Curve{1,2,3,4} = N_W;

Extrude {0,0,H}{
	Surface {1};
	Layers{N_H};
	Recombine;
}

// Surface numbering follows Gmsh's Extrude convention (same pattern as
// tall_cavity.geo): {1}=bottom(hot), {13,17,21,25}=the four side walls in
// extrusion order, {26}=top(cold).
Physical Surface("hot",1)  = {1};
Physical Surface("cold",2) = {26};
Physical Surface("side",3) = {13,17,21,25};
Physical Volume("fluid",4) = {1};
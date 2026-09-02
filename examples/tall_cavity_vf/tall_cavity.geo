H=1.0;
L=H/20.0;

Point(1) = {0,0,0};
Point(2) = {H/2,0,0};
Point(3) = {H/2,H,0};
Point(4) = {0,H,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Transfinite Curve{1,2,3,4} = 21;

Extrude {0,0,L}{
	Surface {1};
	Layers{20};
	Recombine;
}

Physical Surface("hot",1) = {1} ;
Physical Surface("cold",2) = {26};
Physical Surface("topandbot", 3) = {13,21};
Physical Surface("p1", 4) = {17};
Physical Surface("p2", 5) = {25};
Physical Volume("fluid",4) = {1};


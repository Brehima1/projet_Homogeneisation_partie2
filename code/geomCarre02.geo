Mesh.MshFileVersion = 2.2;
// definition du pas du maillage
h = 0.15;
L = 2;
// définition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};
// définition des segments qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// définition des contours fermés
Line Loop(1) = {1,2,3,4};
// définition des surfaces à partir contours fermés
Plane Surface(1) = {1};
// définition des éléments physiques : pour ces éléments, nous pourrons récupérer
//									   les références 
Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {1};

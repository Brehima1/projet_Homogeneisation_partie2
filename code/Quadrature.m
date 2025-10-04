% Quadrature pour approximer A(x, y) sur le triangle (ici, utilisant trois points de quadrature)
function I = Quadrature(F,S1, S2, S3)
% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

%T_chapeau: S1_chapeau, S2_chapeau, S3_chapeau = [0,0] [1,0] [0,1]
%matrice de changement de base
P = [x2-x1, x3-x1; y2-y1, y3-y1];
%vecteur de changement affine
Q = [x1; y1];

% D est, au signe pres, deux fois l'aire du triangle
D = det(P);
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end;

%matrice de changement de base inverse
P_inv = [y3-y1, x1-x3; y1-y2, x2-x1]/D;



% Points de Gauss-Lobatto
M1 = [1/3 1/3];
M2 = [1/5 1/5];
M3 = [1/5 3/5];
M4 = [3/5 1/5];
% Poids de Gauss-Lobatto
w1 = -9/32;
w2 = 25/96;
w3 = 25/96;
w4 = 25/96;
% Calcul de l'integrale
I = abs(D)*(w1*F(P*M1'+Q)+w2*F(P*M2'+Q)+w3*F(P*M3'+Q)+w4*F(P*M4'+Q));
end;


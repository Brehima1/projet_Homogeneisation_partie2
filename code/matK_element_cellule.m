function [Kel] = matK_element_cellule(mat_A_partie2, S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% Calcul de la matrice de raideur élémentaire en P1 Lagrange
%
% SYNOPSIS [Kel] = mat_elem(mat_A_partie2, S1, S2, S3)
%
% INPUT * mat_A_partie2 : fonction qui calcule la matrice A en fonction de x et y
%        S1, S2, S3 : coordonnées des 3 sommets du triangle (vecteurs réels 1x2)
%
% OUTPUT - Kel : matrice de raideur élémentaire (matrice 3x3)
%
% NOTE (1) le calcul est approximé par quadrature à 4 points de Gauss-Lobatto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Préliminaires pour faciliter la lecture :
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Les 3 normales à l'arête opposée (de la longueur de l'arête)
normales = zeros(3, 2);
normales(1, :) = [y2 - y3, x3 - x2];
normales(2, :) = [y3 - y1, x1 - x3];
normales(3, :) = [y1 - y2, x2 - x1];

% D est, au signe près, deux fois l'aire du triangle
D = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
if abs(D) <= eps
  error('L''aire d''un triangle est nulle !');
end

% Calcul de la matrice de raideur
% -------------------------------
cst = 1.0 / (D ^ 2);
Kel = zeros(3, 3);

for i = 1:3
  % Gradient de w_i
  n_i = normales(i, :);

  for j = 1:3
    % Gradient de w_j
    n_j = normales(j, :);
    % Fonction F à intégrer
    F = @(M) cst * (n_i * mat_A_partie2(M(1), M(2)) * n_j');
    % Calcul de l'intégrale avec la méthode de quadrature
    Kel(i, j) = Quadrature(F, S1, S2, S3);
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


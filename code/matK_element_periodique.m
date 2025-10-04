function [Kel] = matK_element_periodique(fonc_A,S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(func_A, S1, S2, S3)
%
% INPUT *  func_A, S1, S2, S3 : la fonction A et les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est aproxime par quadrature a 4 points de Gauss-Lobatto
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
normales = zeros(3, 2);
normales(1, :) = [y2-y3, x3-x2];
normales(2, :) = [y3-y1, x1-x3];
normales(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end;



% calcul de la matrice de raideur
% -------------------------------
cst = 1.0 / (D ^2);
Kel = zeros(3,3);
for i=1:3
  % gradient de w_i
  n_i = normales(i, :);
  for j=1:3
	% A COMPLETER
  n_j = normales(j, :);
  %fonction F a integrer
  F = @(M) cst * (fonc_A(M(1),M(2)) * n_i' )' * n_j';
  Kel(i,j) = Quadrature(F, S1, S2, S3);
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

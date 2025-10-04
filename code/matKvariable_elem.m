function [Kel1] = matKvariable_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par
%          les coordonnees barycentriques
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [-1,-1];
norm(2, :) = [1,0];
norm(3, :) = [0,1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end;

%matrice jacobienne de F
dF = [x2-x1,x3-x1;y2-y1,y3-y1];
idFt = inv(dF');

%integrale de A(x,y)
A_int = -9/32*(dF(1,1)*1/3+dF(1,2)*1/3+x1,dF(2,1)*1/3+dF(2,2)*1/3+y1)+...
        25/96*(A(dF(1,1)*1/5+dF(1,2)*1/5+x1,dF(2,1)*1/5+dF(2,2)*1/5+y1)+...
        A(dF(1,1)*1/5+dF(1,2)*3/5+x1,dF(2,1)*1/5+dF(2,2)*3/5+y1)+...
        A(dF(1,1)*3/5+dF(1,2)*1/5+x1,dF(2,1)*3/5+dF(2,2)*1/5+y1));

% calcul de la matrice de raideur
% -------------------------------
Kel1 = zeros(3,3);
for i=1:3
  for j=1:3
	% A COMPLETER
    Kel1(i,j) = (A_int*(idFt*[norm(i,1); norm(i,2)]))'*(idFt*[norm(j,1); norm(j,2)])*abs(D);
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

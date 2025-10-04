function [erreurL2,erreurH1]=principal1Kvar_neumann(nom_maillage)
% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  %Numeros du triangle
  tri = Numtri(l,:);

  %sommets
  S1=Coorneu(tri(1),:);
  S2=Coorneu(tri(2),:);
  S3=Coorneu(tri(3),:);
  % calcul des matrices elementaires du triangle l
  Kel1=matKvariable_elem(S1, S2, S3);

  Mel=matM_elem(S1, S2, S3);

  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
   for i=1:3
     I = tri(i);
     for j=1:3
       J = tri(j);
       MM(I,J) +=  Mel(i,j);
       KK(I,J) +=  Kel1(i,j);
     endfor
   endfor


end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(: ,1), Coorneu(: ,2)) ;
LL = MM*FF;

% inversion
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));

affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann1 -exacte%s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
erreurL2 = log(sqrt(UU'*MM*UU+UU_exact'*MM*UU_exact-2*UU'*MM*UU_exact)/sqrt(UU_exact'*MM*UU_exact));
D = sqrt(UU'*KK*UU+UU_exact'*KK*UU_exact-2*UU'*KK*UU_exact);
T  = sqrt(UU_exact'*KK*UU_exact);
erreurH1 = log(D/T);
fprintf("erreur L2 = %f et erreur H1 = %f\n", erreurL2, erreurH1);
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%)=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

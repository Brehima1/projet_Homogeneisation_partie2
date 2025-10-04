function [erreurL2,erreurH1]=principal_dirichlet(nom_maillage)

% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
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
  tri = Numtri(l,:);
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
        MM(I,J) += Mel(i,j);
        KK(I,J) += Kel1(i,j);
      end;
  end;
end; % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(: , 1) , Coorneu(: , 2));
LL = MM*FF;

% Projection sur l espace V_0
% ———————————————————
% matrice de projection
num = 1;
PP = zeros(Nbpt-Nbaretes,Nbpt);
for i=1:Nbpt
  if Refneu(i)==0
    PP(num,i) = 1;
    num += 1;
  endif
endfor
AA = MM+KK;
AA0 = PP*AA*PP';
LL0 = PP*LL

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU =PP'*UU0;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('dirichletexact - %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
erreurL2 = log(sqrt(UU'*MM*UU+UU_exact'*MM*UU_exact-2*UU'*MM*UU_exact)/sqrt(UU_exact'*MM*UU_exact));
% Calcul de l erreur H1
erreurH1 = log(sqrt(UU'*KK*UU+UU_exact'*KK*UU_exact-2*UU'*KK*UU_exact)/sqrt(UU_exact'*KK*UU_exact));
fprintf("erreur L2 = %f et erreur H1 = %f\n", erreurL2, erreurH1);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


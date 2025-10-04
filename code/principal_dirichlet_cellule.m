
% =====================================================
%
% Une routine pour la mise en oeuvre des éléments finis P1 Lagrange
% pour l'équation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u = f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================

% Lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre02.msh';
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage);

% ----------------------
% Calcul des matrices EF
% ----------------------

% Déclarations
% ------------
KK = sparse(Nbpt, Nbpt);  % matrice de rigidité
MM = sparse(Nbpt, Nbpt);  % matrice de masse
LL = zeros(Nbpt, 1);      % vecteur second membre

% Boucle sur les triangles
% ------------------------
for l = 1:Nbtri
    % Coordonnées des sommets du triangle
    tri = Numtri(l, :);
    S1 = Coorneu(tri(1), :);
    S2 = Coorneu(tri(2), :);
    S3 = Coorneu(tri(3), :);

    % Calcul des matrices élémentaires du triangle l
    Kel = matK_element_cellule(@mat_A_partie2, S1, S2, S3);
    Mel = matM_elem(S1, S2, S3);

    % Assemblage de la matrice globale et du second membre
    for i = 1:3
      I = tri(i);
      for j = 1:3
        J = tri(j);
        MM(I, J) = MM(I, J) + Mel(i, j);  % Matrice de masse
        KK(I, J) = KK(I, J) + Kel(i, j);  % Matrice de rigidité
      end;
    end;
end ; % fin de la boucle sur les triangles

% Calcul du second membre L
% -------------------------
FF = f_partie2(Coorneu(:, 1), Coorneu(:, 2));  % À adapter en fonction de la définition de f_partie2
LL = MM * FF;

% Projection sur l'espace V_0
% ---------------------------
% Matrice de projection
num = 1;
PP = zeros(Nbpt - Nbaretes, Nbpt);
for i = 1:Nbpt
    if Refneu(i) == 0
        PP(num, i) = 1;
        num = num + 1;
    end;
end;

AA = MM + KK;
AA0 = PP * AA * PP';
LL0 = PP * LL;

% Inversion
% ---------
UU0 = AA0 \ LL0;

% Expression de la solution dans toute la base
% -------------------------------------------
UU = PP' * UU0;

% Visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('dirichlet- %s', nom_maillage));

% Validation
% ----------
validation = 'oui';
if strcmp(validation, 'oui')
    UU_exact = sin(pi * Coorneu(:, 1)) .* sin(pi * Coorneu(:, 2));
    affiche(UU_exact, Numtri, Coorneu, sprintf('DirichletExact - %s', nom_maillage));

    % Calcul de l'erreur L2
    erreurL2 = sqrt((UU - UU_exact)' * MM * (UU - UU_exact) / (UU_exact' * MM * UU_exact));

    % Calcul de l'erreur H1
    erreurH1 = sqrt((UU - UU_exact)' * KK * (UU - UU_exact) / (UU_exact' * KK * UU_exact));

    fprintf("Erreur L2 = %f et erreur H1 = %f\n", erreurL2, erreurH1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        Fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


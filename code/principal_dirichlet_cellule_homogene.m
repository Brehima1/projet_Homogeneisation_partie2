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

u_exact = @(x,y)sin(pi*x).*sin(2*pi*y);
uxx_exact = @(x,y) -pi^2*sin(pi*x).*sin(2*pi*y);
uyy_exact = @(x,y)-4*pi^2*sin(pi*x).*sin(2*pi*y);
uxy_exact = @(x,y) 2*pi^2*cos(pi*x).*cos(2*pi*y);

%Exo2 Question11
%Exo 3 Question 1,2,3,4
A= cell(1, 5);

A{1} = @(x,y) 1;
A{2} = @(x,y) [1, 0; 0, 2];
A{3} = @(x,y)[sqrt(3), 0; 0, 4];
A{4} = @(x,y)[sqrt(3), 0; 0, 4] ;
A{5} = @(x,y)[4*sqrt(3), 0; 0, 2*sqrt(3)]

f= cell(1, 5);

f{1} = @(x,y) 3*pi^2*sin(pi*x).*sin(pi*2*y); % A = I;
f{2} = @(x,y) 9*pi^2*sin(pi*x).*sin(pi*2*y); % A = [1, 0; 0, 2];
f{3} = @(x,y) -2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(2*pi*y) ...
     + (18*pi^2+pi^2*sin(2*pi*x)).*sin(pi*y).*sin(2*pi*y); % A=  [2 + sin(2*pi*x), 0; 0, 4];
f{4} = @(x,y) -2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y)+...
       pi^2*(5*sin(2*pi*x)+ 18).*sin(pi*x).*sin(2*pi*y); % A = [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];
f{5} = @(x,y) @(x,y) -2*pi^2*(sin(2*pi*x)+ 4).*cos(2*pi*x).*cos(pi*y).*sin(2*pi*x)+...
    5*pi^2*(sin(2*pi*x)+ 2).*(sin(2*pi*y)+ 4).*sin(pi*x).*sin(2*pi*y)...
     -4*pi^2*(sin(2*pi*x)+ 2).*sin(pi*x).*cos(2*pi*y).*cos(2*pi*y); % A = (2 + sin(2*pi*x)).*(4 + sin(2*pi*y));


t = 3; % on choisit le probleme a resoudre parmi les 5

validation = 'non';
validation_osc = 'oui'; % Question 1,2,3,4


fonc_A = A{t};
fonc_f = f{t};

% calcul du du tenseur homogénéisé
% --------------------------------
fprintf('Calcul du tenseur homogénéisé\n');
%valeur de eta
eta = 5^(-9);

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per_cellule006.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite

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

   Kel=matK_element_periodique(fonc_A,S1, S2, S3);

   Mel=matM_elem(S1, S2, S3);

  % On fait l'assemblage de la matrice globale et du second membre
  % A COMPLETER
   for i=1:3
      I = tri(i);
      for j=1:3
          J = tri(j);
          MM(I,J) += Mel(i,j);
          KK(I,J) += Kel(i,j);
      end;
  end;

end;% for l

% Calcul du second membre L
%(x1,y1),...,(xn,yn)] est la liste des cornneu, alors XX = [x1,...,xn] et YY = [y1,...,yn]
XX = Coorneu(:,1);
YY = Coorneu(:,2);

%e1 = grad XX et e2 = grad YY, LLX = -K*XX et LLY = -K*YY
LLX = -KK*XX;
LLY = -KK*YY;

% Projection sur l espace V_p
% ———————————————————
% matrice de projection
PP = sparse(Nbpt,Nbpt);
Nbpt_interieur = 0;
for indice = 1:Nbpt
    if Refneu(indice) == 0
        PP(indice,indice) = 1;
        Nbpt_interieur += 1;
    end
end %for indice
Nbpt_bord = Nbpt - Nbpt_interieur;
PP(1,1) = 1;
PP(1,2) = 1;
PP(1,3) = 1;
PP(1,4) = 1;
Nbpt_bord_1 = 0;
for indice = 5:Nbpt_bord
    if Coorneu(indice,2) == 0
       Nbpt_bord_1 += 1;
    else
        break;
    end; %if
end;  %for indice

Nbpt_bord_2 = Nbpt_bord/2 - Nbpt_bord_1 - 2;
for indice = 5: 4 + Nbpt_bord_1
    PP(indice,indice) = 1;
    PP(5 + Nbpt_bord_1 - indice + 4,indice + Nbpt_bord_1 + Nbpt_bord_2) = 1;

end; %for indice

for indice = 5+Nbpt_bord_1:4 + Nbpt_bord_1+Nbpt_bord_2
    PP(indice,indice) = 1;
    PP(5 + Nbpt_bord_1 + Nbpt_bord_2 + Nbpt_bord_1 - indice + 4,indice + Nbpt_bord_1 + Nbpt_bord_2) = 1;
end;


AA = eta*MM+KK;
AAp = PP*AA*PP';
LLXp = PP*LLX;
LLYp = PP*LLY;

% inversion
% ----------
WXp = AAp\LLXp;
WYp = AAp\LLYp;

% Expression de la solution dans toute la base
% ———————
W1 = PP'*WXp;
W2 = PP'*WYp;


% visualisation
% -------------
affiche(W1, Numtri, Coorneu, sprintf('Cellule i=1 - %s', nom_maillage));
affiche(W2, Numtri, Coorneu, sprintf('Cellule i=2 - %s', nom_maillage));

%tenseur homogeneise
Aeff = zeros(2,2);

Aeff(1,1) = (XX + W1)' * (KK * (XX + W1));
Aeff(1,2) = (XX + W1)' * (KK * (YY + W2));
Aeff(2,1) = (YY + W2)' * (KK * (XX + W1));
Aeff(2,2) = (YY + W2)' * (KK * (YY + W2));

Aeff
fprintf('fin du calcul du tenseur homogénéisé\n');
% fin du calcul du tenseur homogénéisé
% ------------

% calcul de la solution homogénéisée
% ----------------------------------
fprintf('Calcul de la solution homogénéisée\n');

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre02.msh';
fprintf('%s\n', nom_maillage);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
KK_std = sparse(Nbpt,Nbpt); % matrice de rigidite standard
fonc_Aeff = @(x,y) Aeff;
fonction_f = @(x,y) -(Aeff(1,1)*uxx_exact(x,y) + ...
                     Aeff(1,2)*uxy_exact(x,y) +...
                     Aeff(2,1)*uxy_exact(x,y) + ...
                     Aeff(2,2)*uyy_exact(x,y));

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
  Kel=matK_element_periodique(fonc_Aeff,S1, S2, S3);
  Kel_std = matK_elem(S1, S2, S3);
  Mel=matM_elem(S1, S2, S3);
  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
  for i=1:3
    I = tri(i);
    for j=1:3
      J = tri(j);
      MM(I,J) += Mel(i,j);
      KK(I,J) += Kel(i,j);
      KK_std(I,J) += Kel_std(i,j);
      end;
  end;

end; % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = fonction_f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% Projection sur l espace V_0
% ———————————————————
% matrice de projection
Point_interieur = find(Refneu == 0);
Nbpt_interieur = size(Point_interieur,1);
PP = sparse(Nbpt_interieur,Nbpt);
for indice = 1:Nbpt_interieur
    PP(indice,Point_interieur(indice)) = 1;
end;
AA = KK;
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UU0;

fprintf('fin du calcul de la solution homogénéisée\n');
% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Dirichlet homogénéisé- %s', nom_maillage));

% validation
% ----------

if strcmp(validation,'oui')
  UU_exact = u_exact(Coorneu(:,1),Coorneu(:,2));
  affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet homogénéisé -exact %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
  errU = UU_exact - UU;
  errU0 = PP*errU;
  errL2 = sqrt(errU' * (MM * errU));
%errL2 = sqrt(errU0' * (PP*MM*PP' * errU0));
% Calcul de l erreur H1
  errH1 = sqrt(errL2^2 + errU' * (KK*errU));
%errH1 = sqrt(errU0' * AA0 * errU0);
% A COMPLETER
  fprintf("erreur L2 = %f et erreur H1 = %f\n", errL2,errH1);
end;

if strcmp(validation_osc,'oui')
  % calcul de la solution ocsillante
  % ---------------------------------
  %question 1,2 ,3 ,4
  eps = 5.^[0:1:7]; % Valeurs de epsilon

  % Initialisation des erreurs
  errL2 = [];
  errH1 = [];

  A_eps = cell(1, numel(eps));
  for h = 1:numel(eps)
    epsilon = eps(h);
    A_eps{h} = @(x, y) fonc_A(x / epsilon, y / epsilon);
  end;
  fprintf('Calcul de la solution oscillante pour différentes valeurs de epsilon\n');
  for h = 1:numel(eps)
    epsilon = eps(h);
    fprintf('epsilon = %.2e\n', epsilon);
    % Initialisation des matrices pour chaque epsilon
    KK = sparse(Nbpt, Nbpt); % Matrice de rigidité
    LL = zeros(Nbpt, 1);     % Second membre
    A_osc = A_eps{h};

    % Assemblage de la matrice KK et calcul du second membre LL
    for l = 1:Nbtri
        tri = Numtri(l, :);
        S1 = Coorneu(tri(1), :);
        S2 = Coorneu(tri(2), :);
        S3 = Coorneu(tri(3), :);
        Kel = matK_element_periodique(A_osc, S1, S2, S3);
        for i = 1:3
            I = tri(i);
            for j = 1:3
                J = tri(j);
                KK(I, J) = KK(I, J) + Kel(i, j);
            end;
        end;
    end;

    % Calcul du second membre FF
    FF = fonction_f(Coorneu(:, 1), Coorneu(:, 2));
    LL = MM * FF;

    % Projection sur l'espace \( V_0 \)
    Point_interieur = find(Refneu == 0);
    Nbpt_interieur = numel(Point_interieur);
    PP = sparse(Nbpt_interieur, Nbpt);
    for indice = 1:Nbpt_interieur
        PP(indice, Point_interieur(indice)) = 1;
    end;

    % Matrice projetée et résolution
    AA = KK;
    AA0 = PP * AA * PP';
    LL0 = PP * LL;
    UU0 = AA0 \ LL0;

    % Retour à l'espace global
    UU_osc = PP' * UU0;

    % Validation des erreurs
    errU_osc = UU - UU_osc;
    errL2_osc = sqrt(errU_osc' * (MM * errU_osc));
    errH1_osc = sqrt(errL2_osc^2 + errU_osc' * (KK_std * errU_osc));
    fprintf("Erreur L2 = %.2e et erreur H1 = %.2e\n", errL2_osc, errH1_osc);

    % Stockage des erreurs
    errL2 = [errL2, errL2_osc];
    errH1 = [errH1, errH1_osc];
end;

% Affichage des erreurs en fonction de epsilon
figure;
hold on;

% Tracé des erreurs en norme L2
log_eps = -log(eps); % Abscisse
log_errL2 = log(errL2); % Ordonnée pour L2
log_errH1 = log(errH1); % Ordonnée pour H1

% Tracé des courbes de l'erreur
plot(log_eps, log_errL2, 'or-', 'LineWidth', 1.5, 'DisplayName', 'Erreur L2');
plot(log_eps, log_errH1, '*b-', 'LineWidth', 1.5, 'DisplayName', 'Erreur H1');

% Ajustements pour le graphique
xlabel('-log(\epsilon)');
ylabel('log(\text{erreur})');
title('Évolution des erreurs L2 et H1 en fonction de -log(\epsilon)');
legend('show', 'Location', 'SouthWest');
grid on;
hold off;
end;

% =====================================================
%
% Routine pour les EF P1 Lagrange
% Résolution de -div(A grad u) + u = f avec conditions périodiques
% sur un maillage donné
%
% =====================================================

clear; clc;

% Fonctions exactes pour validation
func_w1_exact = @(x, y) zeros(size(x));
func_w2_exact = @(x, y) zeros(size(x));

% Validation activée/désactivée
validation = 'oui'; % 'oui' pour activer la validation

% Paramètre eta
eta = 2^(-8);

% Liste des fonctions A selon les exercices
if strcmp(validation, 'oui')
    % Cas pour validation
    func_A_list = {
        @(x, y) 1;
        @(x, y) [1, 0; 0, 2]
    };
else
    % Cas général
    func_A_list = {
        @(x, y) 1;
        @(x, y) [1, 0; 0, 2];
        @(x, y) [2 + sin(2 * pi * x), 0; 0, 4];
        @(x, y) [2 + sin(2 * pi * x), 0; 0, 4 + sin(2 * pi * x)];
        @(x, y) (2 + sin(2 * pi * x)) .* (4 + sin(2 * pi * y))
    };
end

% Chargement du maillage
nom_maillage = 'geomCarre_per0075.msh';
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage);

% Boucle sur les différents cas de A
for index_func = 1:length(func_A_list)
    func_A = func_A_list{index_func};

    % Initialisation des matrices globales
    KK = sparse(Nbpt, Nbpt); % Matrice de rigidité
    MM = sparse(Nbpt, Nbpt); % Matrice de masse

    % Assemblage des matrices EF
    for l = 1:Nbtri
        % Sommets du triangle courant
        tri = Numtri(l, :);
        S1 = Coorneu(tri(1), :);
        S2 = Coorneu(tri(2), :);
        S3 = Coorneu(tri(3), :);

        % Calcul des matrices élémentaires
        Kel = matK_elem_dynamique(func_A, S1, S2, S3);
        Mel = matM_elem(S1, S2, S3);

        % Assemblage dans les matrices globales
        for i = 1:3
            I = tri(i);
            for j = 1:3
                J = tri(j);
                MM(I, J) = MM(I, J) + Mel(i, j);
                KK(I, J) = KK(I, J) + Kel(i, j);
            end
        end
    end

    % Second membre
    XX = Coorneu(:, 1);
    YY = Coorneu(:, 2);
    LLX = -KK * XX;
    LLY = -KK * YY;

    % Construction de la matrice de projection
    PP = sparse(Nbpt, Nbpt);
    Nbpt_interieur = sum(Refneu == 0);

    % Points intérieurs
    PP(Refneu == 0, Refneu == 0) = eye(Nbpt_interieur);

    % Conditions périodiques
    Nbpt_bord = Nbpt - Nbpt_interieur;
    for i = 1:Nbpt_bord / 2
        PP(Nbpt_interieur + i, Nbpt_interieur + Nbpt_bord / 2 + i) = 1;
    end

    % Assemblage matrice et second membre projetés
    AA = eta * MM + KK;
    AAp = PP * AA * PP';
    LLXp = PP * LLX;
    LLYp = PP * LLY;

    % Résolution du système
    WXp = AAp \ LLXp;
    WYp = AAp \ LLYp;

    % Expression dans la base complète
    W1 = PP' * WXp;
    W2 = PP' * WYp;

    % Visualisation
    affiche(W1, Numtri, Coorneu, sprintf('Cellule i=1 - %s', nom_maillage));
    affiche(W2, Numtri, Coorneu, sprintf('Cellule i=2 - %s', nom_maillage));

    % Calcul du tenseur homogénéisé
    Aeff = zeros(2, 2);
    Aeff(1, 1) = (XX + W1)' * (KK * (XX + W1));
    Aeff(1, 2) = (XX + W1)' * (KK * (YY + W2));
    Aeff(2, 1) = (YY + W2)' * (KK * (XX + W1));
    Aeff(2, 2) = (YY + W2)' * (KK * (YY + W2));

    fprintf('Cas %d\n', index_func - 1);
    disp(Aeff);

    % Validation (si activée)
    if strcmp(validation, 'oui')
        W1_exact = func_w1_exact(Coorneu(:, 1), Coorneu(:, 2));
        W2_exact = func_w2_exact(Coorneu(:, 1), Coorneu(:, 2));

        % Calcul des erreurs
        errX = W1 - W1_exact;
        errY = W2 - W2_exact;
        errL2 = sqrt(errX' * (MM * errX) + errY' * (MM * errY));
        errH1 = sqrt(errL2^2 + errX' * (KK * errX) + errY' * (KK * errY));

        fprintf("Erreur L2 = %f et Erreur H1 = %f\n", errL2, errH1);
    end
end


function val = mat_A_partie2(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%Partie2  du TP
eps =100;

%val = [1 0;0 1];
%val = [1, 0; 0, 2];% i)question4
%val =  [2 + sin(2*pi*x/eps), 0; 0, 4]; %ii)question4
%val = [2 + sin(2*pi*x/eps), 0; 0, 4 + sin(2*pi*y/eps)]; % iii)question4
val = (2 + sin(2*pi*x/eps)).*(4 + sin(2*pi*y/eps))*[1 0;0 1]; % iv)question4





end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

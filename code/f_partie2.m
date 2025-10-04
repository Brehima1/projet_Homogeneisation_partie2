function val = f_partie2(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%Partie2 du TP
%val = 2*pi^2*sin(pi*x).*sin(pi*y); % A = I
%val = 3*pi^2*sin(pi*x).*sin(pi*y); % A = [1, 0; 0, 2]; i)question4
%val =-2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y) ...
     %+ pi^2*(sin(2*pi*x)+ 2).*sin(pi*x).*sin(pi*y)...
     %+ 4*pi^2*sin(pi*x).*sin(pi*y); % A=  [2 + sin(2*pi*x), 0; 0, 4]; ii)question4
%val = -2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y) ...
      %+ pi^2*(sin(2*pi*x)+ 4).*sin(pi*x).*sin(pi*y) ...
      %+ pi^2*(sin(2*pi*x)+ 2).*sin(pi*x).*sin(pi*y); % A = [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];  iii)question4
val =-2*pi^2*(sin(2*pi*x)+ 2).*cos(2*pi*y).*cos(pi*y).*sin(pi*x) ...
     - 2*pi^2*(sin(2*pi*y)+ 4).*cos(2*pi*x).*cos(pi*x).*sin(pi*y)...
    + 2*pi^2*(sin(2*pi*x)+ 2).*(sin(2*pi*y)+ 4).*sin(pi*x).*sin(pi*y); % A = (2 + sin(2*pi*x)).*(4 + sin(2*pi*y))*Id;  iv)question4

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

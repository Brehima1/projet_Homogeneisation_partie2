%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecture_msh.m:
% routine de lecture de fichiers de maillages triangulaires 2D au format .msh
%
% SYNOPSIS [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nomfile)
%
% INPUT * nom_maillage : le nom d'un fichier de maillage au format msh
%                   AVEC SON SUFFIXE .msh
%
% OUTPUT - Nbpt : nbre de sommets (entier)
%        - Nbtri : nbre de triangles (entier)
%        - Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        - Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        - Numtri : liste de triangles
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        - Reftri : reference des triangles (matrice entiere Nbtri x 1)
%        - Nbaretes : Nombre d'aretes sur la frontiere.
%        - Numaretes(Nbaretes,2) = numero des 2 Noeuds de chaque Arete de la frontiere
%		 - Refaretes(Nbaretes,1) = Reference de chaque ar?te de la
%		 fronti?re
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020
function [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nomfile)
fid=fopen(nomfile,'r');
if fid <=0
   msg=['Le fichier de maillage : ' nomfile ' n''a pas ?t? trouv?'];
   error(msg);
end
while ~strcmp(fgetl(fid),'$Nodes'), end
Nbpt = str2num(fgetl(fid));
Coorneu = zeros(Nbpt,2);
Refneu = zeros(Nbpt,1);
Nbaretes = 0;
Numaretes = [];
Refaretes = [];
RefneuBis = zeros(Nbpt,1);
for i=1:Nbpt
tmp= str2num(fgetl(fid));
Coorneu(i,:) = tmp(2:3);
end
while ~strcmp(fgetl(fid),'$Elements'), end
Nbtri = str2num(fgetl(fid));
tmp= str2num(fgetl(fid));test = tmp(2);
while ( (test~=15) && (test~=1) && (test~=2)), tmp= str2num(fgetl(fid)); test = tmp(2); Nbtri = Nbtri-1;end
% on traite les noeuds des coins
while ((test~=1) && (test~=2))
     Refneu(tmp(end))=tmp(end-2);
     tmp= str2num(fgetl(fid)); test = tmp(2); Nbtri = Nbtri-1;
 end
% on traite les noeuds du bord
k=1;
while (test~=2)
    RefneuBis(tmp(6:7)')=tmp(4);
    Numaretes= [Numaretes;tmp(6:7)];
    Refaretes= [Refaretes;tmp(4)];
    tmp= str2num(fgetl(fid));
    test = tmp(2); Nbaretes=Nbaretes+1; Nbtri=Nbtri-1;k=k+1;
end
Refneu(find(Refneu==0))=RefneuBis(find(Refneu==0));
Numtri = zeros(Nbtri,3);
Reftri = zeros(Nbtri,1);
% on traite les references des triangles
for i=1:Nbtri
    Numtri(i,:) = tmp(end-2:end);
    Reftri(i)=tmp(end-3);
    tmp= str2num(fgetl(fid));
end
fclose(fid);

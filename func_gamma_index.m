function [g_index,pass]=func_gamma_index(z_measure,dose_measure,z_simulation,dose_simulation,dta,dose_diff)

R= [z_simulation dose_simulation];   %calcul;
E= [z_measure dose_measure];   %mesure;

G=zeros(size(z_simulation));
k=0;
t=0;
maxMatrice=max(z_simulation);
%============================== 1
% Mise à l'echelle
%==============================
Mmaxi=max(dose_measure);
Cmaxi=max(dose_simulation);

%z_simulation=transpose(z_simulation);  dose_simulation=transpose(dose_simulation);
%Nmes=dose_measure*100/Mmaxi; 
Nmes=dose_measure*100/Mmaxi;
%Ncal=dose_simulation*100/Mmaxi;
%Ncal=dose_simulation*100/Cmaxi;
Ncal=dose_simulation*100/Cmaxi;

R=[z_measure Nmes];
E=[z_simulation Ncal];


%==============================+
% calcul indice gamma
%==============================
DTA=dta; M=dose_diff;
%if length(z_simulation)>length(z_measure)
%    run=length(z_measure);
%
    run=length(z_simulation);
%end
%display('size(z_simulation)');    display(size(z_simulation));
%display('size(z_measure)');    display(size(z_measure));
for i=1+2:run-2;
   
       t=t+1;
        V1=sqrt((1/DTA)^2+((E(i,2)-R(i+1,2))/M)^2);
        V2=sqrt((1/DTA)^2+((E(i,2)-R(i-1,2))/M)^2);
        V3=sqrt((0/DTA)^2+((E(i,2)-R(i,2))/M)^2);
        V4=sqrt((2/DTA)^2+((E(i,2)-R(i-2,2))/M)^2);
        V5=sqrt((2/DTA)^2+((E(i,2)-R(i+2,2))/M)^2);
        
     V=[V1,V2,V3,V4,V5];
        G(i)=min(V) ;
        
       if G(i)<1;
        k=k+1;
       end
        
end
%display(G);

g_index=G;
passage=k/t*100;
pass=passage;
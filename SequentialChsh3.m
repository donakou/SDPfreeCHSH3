clear



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  suivre sur le tableau excel pour mieux comprendre%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%matrice définissant le polynôme

M1BA= zeros(49,49);

M1BA(2,14)=1; M1BA(2,32)=1; M1BA(2,20)=-1; M1BA(2,44)=-1;

M1BA(3,21)=1; M1BA(3,39)=1; M1BA(3,27)=-1; M1BA(3,33)=-1;

M1BA(4,28)=1; M1BA(4,46)=1; M1BA(4,16)=-1; M1BA(4,40)=-1;

M1BA(5,23)=1; M1BA(5,35)=1; M1BA(5,17)=-1; M1BA(5,41)=-1;

M1BA(6,42)=1; M1BA(6,30)=1; M1BA(6,24)=-1; M1BA(6,48)=-1;

M1BA(7,19)=1; M1BA(7,49)=1; M1BA(7,31)=-1; M1BA(7,37)=-1;


X=sdpvar(49,49,'hermitian')

constraints= [X>=0];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Ancienne matrice des moments d'ordre%1 et 2
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Les mesures et leurs projecteurs

mesure= [[2,3,4];[5,6,7];[8,9,10];[11,12,13]];

       


% coefficients à 0

for k = [1,2,3,4]  % pour chacune des mesures
    for i= mesure(k,:)
        for j= mesure(k,:)
            if not(i==j)
               constraints=[constraints, X(i,j)==0 ];
            end
        end
    end
end
  


% coefficients à 1
 constraints=[constraints,X(1,1)==1];


 
 % contrainte d'égalité des coefficients
 for i=2:13
     constraints=[constraints,X(1,i) == X(i,i) ];
 end



%somme des probas presqu'à 1 

for k = [1,2,3,4] , % pour chacune des mesures
    
   constraints=[constraints, X(1,mesure(k,1))+X(1,mesure(k,2))+X(1,mesure(k,3))<=1];
    
end
%le <= pour tenir compte de la dimension plus grande que 3


%autres sommes des probas  

for i = 5:13
  
     constraints=[constraints, X(2,i)+X(3,i)+X(4,i) <= X(1,i)];
end
        
for i = 8:13
  
     constraints=[constraints, X(5,i)+X(6,i)+X(7,i) <= X(1,i)];
end
for i = 1:13
  
     constraints=[constraints, X(8,i)+X(9,i)+X(10,i) <= X(1,i)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Les moments d'ordre 2 apparaissant dans la nouvelle matrice des momments sont 
%égaux à ceux de la matrice précédente

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cpt1=1;  % compteur des projo
cpt2=0;  % compteur des mesures

while (13 + 6*cpt2+ cpt1 <50) % par rapport à la taille de la matrice
    
    for cpt1 =1:6
        %display([1,13 + 6*cpt2+ cpt1,8+cpt2,cpt1+1])
        
        constraints=[constraints,X(1,13 + 6*cpt2+ cpt1)== X(8+cpt2,cpt1+1)];
        
    end
    cpt2=cpt2+1;
    cpt1=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Les moments d'ordre%3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%contraintes sur les sommes


for i= 14 : 49  % Pour chaque colonne
    
    for cpt2 = 0 : 3  % pour chacune des mesures
        %display([i,cpt2*3+2,i,cpt2*3+3,i,cpt2*3+4])
        constraints=[constraints, X(cpt2*3+2,i)+ X(cpt2*3+3,i) + X(cpt2*3+4,i) <= X(1,i)];
        
    end
end



%coefficients nuls

for i = 14 : 19
    
     constraints=[constraints,X(9,i)==0,X(10,i)==0];
end

for i= 20:25
     constraints=[constraints, X(8,i)==0,X(10,i)==0];
end

for i= 26:31
     constraints=[constraints, X(8,i)==0,X(9,i)==0];
end



for i = 32 : 37
    
     constraints=[constraints,X(12,i)==0,X(13,i)==0];
end

for i= 38:43
     constraints=[constraints, X(11,i)==0,X(13,i)==0];
end

for i= 44:49
     constraints=[constraints, X(11,i)==0,X(12,i)==0];
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Les moments d'ordre% 4 sus 3 
                        %(conditions diagonales)
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cpt1=0;

while 14 + cpt1 < 50  %Taille de la matrice

    %on se place sur le 1er coeff qui interesse et on parcoure la matrice
for i=0:5
    for j=0:5
        constraints=[constraints, X(14+cpt1+i,14+cpt1+j) == X(2+i,14+cpt1+j)];
        %display ([14+cpt1+i,14+cpt1+j,2+i,14+cpt1+j])
    end
end

cpt1=cpt1+6 ;
% un peu de redondance, mais plus facile à comprendre
end





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %Les moments d'ordre% 4 
%                         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%         %%%%%%%%%% Relation avec les moments d'ordre 3%%%%%%%%
%         
% for i = 32:49   % sur chacune des colonnes ne commençant pas par des 0
%     for cpt = 0 : 2
%          constraints=[constraints, X(14+6*cpt,i)+X(14+6*cpt+1,i)+X(14+6*cpt+2,i)<= X(8+cpt,i)];
%          constraints=[constraints, X(14+6*cpt+3,i)+X(14+6*cpt+4,i)+X(14+6*cpt+5,i)<= X(8+cpt,i)];
%          
%          %display ([14+6*cpt,i,14+6*cpt+1,i,14+6*cpt+2,i,8+cpt,i ])
%          %display ([14+6*cpt+3,i,14+6*cpt+4,i,14+6*cpt+5,i,8+cpt,i ])
%          
%     end
%     
% end
% 
% 
% 
%         %%%%%%%% Coefficients nuls %%%%%%
%         
% for i = 20 : 31 
%     for j = 14 : 19
%         
%         constraints=[constraints, X(i,j)==0];
%         
%     end
%     
% end
% 
% 
% for i = 26 : 31 
%     for j =20 : 25
%         
%         constraints=[constraints, X(i,j)==0];
%         
%     end
%     
% end
% 
% 
% for i = 38 : 49 
%     for j = 32 : 37
%         
%         constraints=[constraints, X(i,j)==0];
%         
%     end
%     
% end
% 
% for i = 44 : 49 
%     for j = 38 : 43
%         
%         constraints=[constraints, X(i,j)==0];
%         
%     end
%     
% end
% 

%


  
ops = sdpsettings('solver','sedumi')
solvesdp(constraints,-trace(M1BA*X),ops)
double(X)
-(-trace(M1BA*double(X)))

G=double(X)
GG = G(1:13,1:13)
[P,D]=eig(GG)
B=real(P*sqrt(D))



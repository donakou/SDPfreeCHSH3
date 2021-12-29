% Pour Certifier le protocole Gabriel, nous faisons des euristiques 
% permetant de voir à quel point l'on s'éloigne de la stratégie optimale 
% pour l'entropie (c'est à dire des probas pour chaque mesure de 1/3
 

%strategie(1,2,(2/3)*(3+sqrt(3)))

i=1
for b= 0:+0.1:3.1

    H(:,:,i)= securite(b)
    i=i+1
end

%%
function stra = strategie(u,v,bellValue) %u et v sont les positions de la première proba de la mesure
ops = sdpsettings('solver','mosek')

%matrice définissant le polynôme
M1= zeros(13,13)
M1(2,8 : 13) = [1 -1 0 1 0 -1]
M1(3,8 : 13) = [0 1 -1 -1 1 0]
M1(4,8 : 13) = [-1 0 1 0 -1 1]
M1(5,8 : 13) = [-1 1 0 1 -1 0]
M1(6,8 : 13) = [0 -1 1 0 1 -1]
M1(7,8 : 13) = [1 0 -1 -1 0 1]

X=sdpvar(13,13)
solvesdp([X>=0, 
    % coefficients à 0
X(2,3)==0, X(2,4)==0, X(3,4)==0,
X(5,6)==0, X(5,7)==0, X(6,7)==0,
X(8,9)==0, X(8,10)==0, X(9,10)==0, 
X(11,12)==0, X(11,13)==0, X(12,13)==0,

% coefficients à 1
X(1,1)==1,

% contrainte d'égalité des coefficients
X(1,2)==X(2,2) , X(1,3)==X(3,3) , X(1,4)== X(4,4) , X(1,5) == X(5,5) ,
X(1,6) == X(6,6) , X(1,7)==X(7,7), X(1,8)==X(8,8), X(1,9)==X(9,9),
X(1,10)==X(10,10), X(1,11)==X(11,11), X(1,12)==X(12,12), X(1,13)==X(13,13),


%somme des probas à 1 
X(1,2) + X(1,3) + X(1,4)  ==1, X(1,5) + X(1,6) + X(1,7)  ==1,
X(1,8) + X(1,9) + X(1,10) ==1, X(1,11) + X(1,12) + X(1,13)==1 ,

% Remettre tous les coefficients de la matrice positifs car ils sont des probas
%On fait que le triangle supérieur car la matrice est symétrique

X(1,2)>=0, X(1,3)>=0, X(1,4)>=0, X(1,5)>=0, X(1,6)>=0, X(1,7)>=0, X(1,8)>=0, X(1,9)>=0,
X(1,10)>=0, X(1,11)>=0, X(1,12)>=0, X(1,13)>=0,
X(2,3)>=0, X(2,4)>=0, X(2,5)>=0, X(2,6)>=0, X(2,7)>=0, X(2,8)>=0, X(2,9)>=0,
X(2,10)>=0, X(2,11)>=0, X(2,12)>=0, X(2,13)>=0,
X(3,4)>=0, X(3,5)>=0, X(3,6)>=0, X(3,7)>=0, X(3,8)>=0, X(3,9)>=0,
X(3,10)>=0, X(3,11)>=0, X(3,12)>=0, X(3,13)>=0,
X(4,5)>=0, X(4,6)>=0, X(4,7)>=0, X(4,8)>=0, X(4,9)>=0, X(4,10)>=0, X(4,11)>=0, X(4,12)>=0, 
X(4,13)>=0,
X(5,6)>=0, X(5,7)>=0, X(5,8)>=0, X(5,9)>=0, X(5,10)>=0, X(5,11)>=0, X(5,12)>=0, X(5,13)>=0,
X(6,7)>=0, X(6,8)>=0, X(6,9)>=0, X(6,10)>=0, X(6,11)>=0, X(6,12)>=0, X(6,13)>=0,
X(7,8)>=0, X(7,9)>=0, X(7,10)>=0, X(7,11)>=0, X(7,12)>=0, X(7,13)>=0,
X(8,9)>=0, X(8,10)>=0, X(8,11)>=0, X(8,12)>=0, X(8,13)>=0,
X(9,10)>=0, X(9,11)>=0, X(9,12)>=0, X(9,13)>=0,
X(10,11)>=0, X(10,12)>=0, X(10,13)>=0,
X(11,12)>=0, X(11,13)>=0, X(12,13)>=0,

trace(M1*X)== bellValue 

], -X(u,v),ops) % Cela nous donne le cas pire (en terme de probas) 
% permettant d'obtenir la bellValue initialement fixée. 
double(X)
stra = X(u,v)
end


%on repète cette procedure sur les coefficcients de Proba conjointes qui 
%utilisées pour générer les trits aléatoires
function TableauProba = securite(bellValue) 

TableauProba = zeros(13,13)
for k = 2:4
    for l = 5:7
        TableauProba(k,l)=strategie(k,l,bellValue)
        TableauProba(l,k)=TableauProba(k,l) %la matrice est symétrique
    end
end


for k = 8:10
    for l = 11:13
        TableauProba(k,l)=strategie(k,l,bellValue)
        TableauProba(l,k)=TableauProba(k,l) %la matrice est symétrique
    end
end

end

%T=-(-trace(M1*double(X)))
% car sdp(f(x)) donne le min de f(x) et on cherche le max avec la relation
% max(f(x)) = - min(-f(x))

%format longE 
%T

%G=double(X)
%[P,D]=eig(G)
%B=real(P*sqrt(D))
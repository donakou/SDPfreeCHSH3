%matrice symétrique définissant le polynôme
M1= zeros(13,13)
M1(2,8 : 13) = [1 -1 0 1 0 -1]
M1(3,8 : 13) = [0 1 -1 -1 1 0]
M1(4,8 : 13) = [-1 0 1 0 -1 1]
M1(5,8 : 13) = [-1 1 0 1 -1 0]
M1(6,8 : 13) = [0 -1 1 0 1 -1]
M1(7,8 : 13) = [1 0 -1 -1 0 1]

M1(8 : 13,2) = [1 -1 0 1 0 -1]
M1(8 : 13,3) = [0 1 -1 -1 1 0]
M1(8 : 13,4) = [-1 0 1 0 -1 1]
M1(8 : 13,5) = [-1 1 0 1 -1 0]
M1(8 : 13,6) = [0 -1 1 0 1 -1]
M1(8 : 13,7) = [1 0 -1 -1 0 1]

M1=1/2* M1


% coefficients à 0
u23= zeros(13,13)
u23(2,3)= 1

u24= zeros(13,13)
u24(2,4) = 1

u32= zeros(13,13)
u32(3,2)= 1

u34= zeros(13,13)
u34(3,4)= 1

u42= zeros(13,13)
u42(4,2)= 1

u43=zeros(13,13)
u43(4,3)= 1

u56= zeros(13,13)
u56(5,6) = 1

u57 = zeros(13,13)
u57(5,7)= 1

u65 = zeros(13,13)
u65(6,5)=1

u67 = zeros(13,13)
u67(6,7)=1

u75 = zeros(13,13)
u75(7,5)=1

u76 = zeros(13,13)
u76(7,6)=1

u89= zeros(13,13)
u89(8,9)= 1

u810= zeros(13,13)
u810(8,10) = 1

u98= zeros(13,13)
u98(9,8)= 1

u910= zeros(13,13)
u910(9,10)= 1

u108= zeros(13,13)
u108(10,8)= 1

u109=zeros(13,13)
u109(10,9)= 1

u1112= zeros(13,13)
u1112(11,12) = 1

u1113 = zeros(13,13)
u1113(11,13)= 1

u1211 = zeros(13,13)
u1211(12,11)=1

u1213 = zeros(13,13)
u1213(12,13)=1

u1311 = zeros(13,13)
u1311(13,11)=1

u1312 = zeros(13,13)
u1312(13,12)=1



% Matrices des contraintes d'égalité des coefficients
eq0=zeros(13,13)
eq0(1,1)=1


eq1= zeros(13,13)
eq1(1,2)=1
eq1(2,2)=-1

eq2 = zeros(13,13)
eq2(3,1)= 1
eq2(3,3)=-1


eq3= zeros(13,13)
eq3(4,1)=1
eq3(4,4)=-1

eq4=zeros(13,13)
eq4(5,1)=1
eq4(5,5) = -1

eq5=zeros(13,13)
eq5(6,1)=1
eq5(6,6)=-1

eq6=zeros(13,13)
eq6(7,1)=1
eq6(7,7)=-1

eq7 = zeros(13,13)
eq7(8,1)= 1
eq7(8,8)=-1

eq8 = zeros(13,13)
eq8(9,1)= 1
eq8(9,9)=-1


eq9 = zeros(13,13)
eq9(10,1)= 1
eq9(10,10)=-1

eq10 = zeros(13,13)
eq10(1,11)= 1
eq10(11,11)=-1

eq11 = zeros(13,13)
eq11(1,12)= 1
eq11(12,12)=-1

eq12 = zeros(13,13)
eq12(1,13)= 1
eq12(13,13)=-1



% matrices des contraintes des probas et positivité

s1=zeros(13,13)
s1(1,2:4)=[1,1,1]

s2=zeros(13,13)
s2(1,5:7)=[1,1,1]

s3=zeros(13,13)
s3(1,8:10)=[1,1,1]

s4=zeros(13,13)
s4(1,11:13)=[1,1,1]

% On va maintenir ces coeff en positif
% car ils apparaissent négatifs et ce sont des probas
pos1=zeros(13,13)
pos1(2,9)=1

pos2=zeros(13,13)
pos2(2,13)=1

pos3=zeros(13,13)
pos3(3,10)=1

pos4=zeros(13,13)
pos4(3,11)=1

pos5=zeros(13,13)
pos5(4,8)= 1

pos6=zeros(13,13)
pos6(4,12)=1

pos7=zeros(13,13)
pos7(5,8)=1

pos8=zeros(13,13)
pos8(5,12)=1

pos9=zeros(13,13)
pos9(6,9)=1

pos10=zeros(13,13)
pos10(6,13)=1

pos11=zeros(13,13)
pos11(7,10)=1

pos12=zeros(13,13)
pos12(7,11)=1



X=sdpvar(13,13)
solvesdp([X>=0, 
    % coefficients à 0
trace(u23*X)==0, trace(u24*X)==0, 
trace(u32*X)==0,trace(u34*X)==0,
trace(u42*X)==0, trace(u43*X)==0,
trace(u56*X)==0, trace(u57*X)==0, 
trace(u65*X)==0, trace(u67*X)==0,
trace(u75*X)==0, trace(u76*X)==0,
trace(u89*X)==0, trace(u810*X)==0, 
trace(u98*X)==0,trace(u910*X)==0,
trace(u108*X)==0, trace(u109*X)==0,
trace(u1112*X)==0, trace(u1113*X)==0, 
trace(u1211*X)==0, trace(u1213*X)==0,
trace(u1311*X)==0, trace(u1312*X)==0,

% coefficients à 1
trace(eq0*X)-1==0,

% contrainte d'égalité des coefficients
trace(eq1*X)==0,
trace(eq2*X)==0,
trace(eq3*X)==0,
trace(eq4*X)==0,
trace(eq5*X)==0,
trace(eq6*X)==0,
trace(eq7*X)==0,
trace(eq8*X)==0,
trace(eq9*X)==0,
trace(eq10*X)==0,
trace(eq11*X)==0,
trace(eq12*X)==0,

% egalité avec somme des coeffs
%trace(eq13*X)==0,
%trace(eq14*X)==0,
%trace(eq15*X)==0,
%somme des probas
trace(s1*X)==1,
trace(s2*X)==1,
trace(s3*X)==1,
trace(s4*X)==1,
%trace(pos1*X)>=0,
%trace(pos2*X)>=0,
%trace(pos3*X)>=0,
%trace(pos4*X)>=0,
%trace(pos5*X)>=0,
%trace(pos6*X)>=0,
%trace(pos7*X)>=0,
%trace(pos8*X)>=0,
%trace(pos9*X)>=0,
%trace(pos10*X)>=0,
%trace(pos11*X)>=0,
%trace(pos12*X)>=0,
%trace(s5*X)==1,
%trace(s6*X)==1,
],-trace(M1*X))
double(X)
-(-trace(M1*double(X)))
%T=-(-trace(M1*double(X)))
% car sdp(f(x)) donne le min de f(x) et on cherche le max avec la relation
% max(f(x)) = - min(-f(x))

%format longE 
%T

G=double(X)
[P,D]=eig(G)
B=real(P*sqrt(D))

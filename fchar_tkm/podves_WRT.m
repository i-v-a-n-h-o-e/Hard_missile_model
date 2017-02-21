% Функция: матрица B по углам ГСП
function M=podves_WRT(Angles);
% Angles=[fi psi tet]
C=cos(Angles);      S=sin(Angles);
Sfi=[1 0 0;     0  C(1)  S(1);    0  -S(1)  C(1)];
Spsi=[C(2)  0 -S(2);     0 1 0;     S(2)  0  C(2)];
Stet=[C(3)  S(3)  0;  -S(3)  C(3)  0;   0 0 1]; 
M=[0 -1 0; 1 0 0; 0 0 1]*Sfi'*Spsi'*Stet';
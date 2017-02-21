function Om=Om_WRT(psi,tet,fie_t,psi_t,tet_t);
%Функция угл. скор. по производным углов
cp=cos(psi);    sp=sin(psi);
ct=cos(tet);    st=sin(tet);
M=[cp*ct  st  0; -cp*st  ct  0; sp  0  1];
Om=M*[fie_t;psi_t;tet_t];
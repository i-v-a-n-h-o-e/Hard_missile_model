function boderus_A_v2(Name,gr1,gr2)
%Расчет АФЧХ для непрерывной или дискретной LTI-системы
dls = 3;         % Ном. расстояние между точками полярн. графика, дб
nidA = 30;       % Уровень нуля амплитуд полярн. графика, дб
Rzon = 6;        % Радиус окрестности крит. точки
Vnorm = [1 1.25 1.6 2 2.5 3 4 5 6 8 10];  % Ряд чисел для границ шкалы графиков
Nd = 1000;%100;%500;%1000;%40;         % Мин. число точек на декаду частоты
rab = 5^0.5/2;
arrow = 2*[-rab+1i/2;0;-rab-1i/2;-rab+1i/2];
Rez = [];
%--------------------------------------------------------------------------%
Tm = Name.Ts;

if ~isempty(Name.notes)
   tit = Name.notes;
else
   tit = '';
end

%------------Определение omn,omk,om_mx,PS,Nppo,Nppz
P = pole(Name);
PS = length(P);

if Tm > 0   % Случай с дискретной системой
    Pz = pole(absorbDelay(feedback(Name,1,1)));
    PZ = length(Pz);
    D1 = sort(abs(P - 1));
    q = 1; d1 = 0;

    while d1 < (1e-5)
        d1 = D1(q);
        q = q + 1;
    end

    omn = 0.001*d1/Tm;      % Стартовая частота годографа
    omk = pi/Tm;            % Конечная частота годографа
    oma = angle(P)/Tm;

    for s = PS:-1:1
        if abs(P(s)) < 2/3        %2/3!!!!!!!!!!!!!!!!!!!!!!!
            oma(s) = [];
        end
    end

    for s = PZ:-1:1
        if abs(Pz(s)) <= 1
            Pz(s) = [];
        end
    end
    if ~isempty(Pz)
        Tunst = [angle(Pz)/Tm, Tm./log(abs(Pz))];
    else
        Tunst=[];
    end
    Nppo = sum(abs(P)>1);
else    % Tm = 0 - Случай с непрерывной системой
    Pz = pole(pade(feedback(Name,1,1)));
    PZ = length(Pz);
    D1 = sort(abs(P));
    q = 1; d1 = 0;
    while d1 < (1e-5)
        d1 = D1(q);
        q = q + 1;
    end
    omn = 0.001*d1;       %0.6!!!!
    Ima = sort(imag(P));
    omk = Ima(PS)/0.001;       %0.6!!!!
    oma = imag(P);
    for s=PS:-1:1
        if real(P(s)) < -10        %-10!!!!!!!!!!!!!!!!!!!!!!!
            oma(s) = [];
        end
    end
    for s=PZ:-1:1
        if real(Pz(s)) <= 0
             Pz(s) = [];
        end
    end
    if ~isempty(Pz)
        Tunst = [imag(Pz), 1./real(Pz)];
    else
        Tunst = [];
    end
   Nppo = sum(real(P) > 0);
end    %Tm=0

rab1 = sort(oma); 
while rab1(1) < omn
    rab1(1) = [];
end
for s = length(rab1):-1:2
    if (rab1(s)/rab1(s-1) < 1.005)     %1.005!!!!!!!!!!!!!!!!!!!!!!1
        rab1(s) = [];
    end
end
om_mx = [rab1; omk];
Nmx = zeros(size(om_mx));
disp('nmx')
disp(length(om_mx))

%--------------------------------------------------------------------------%
Knd = 10^(1/Nd);  gam2 = 5/57.3; gam1 = gam2/2.5;
Kuv = 1; ome = omn;  n = 1;  gm = 1;

while gm <= length(om_mx)
    if n <= 2
        if n == 2
            ome = ome*Knd;
        end
        G.ome(n) = ome;            % omega, Gout(n,2);
        if Tm > 0
            G1 = evalfr(Name, exp(1i*ome*Tm)); %feval(ProbA,p);
        else
            G1 = evalfr(Name, 1i*ome);
        end
        G.mdb(n) = 20*log10(abs(G1));   %Gout(n,3)
        G.ang(n) = angle(G1);   %Gout(n,4)
        n = n + 1;
    else    %n > 2
        if  G.mdb(n-1) > -28          %-28 
            Ks = G.ome(n-1)/G.ome(n-2);
            omea = Ks^Kuv*G.ome(n-1);
        else
            omea = Knd*G.ome(n-1);
        end
        if ~isempty(Rez)&&(omea >= Rez(1,1))&&(Rez(1,1) > 1.005*G.ome(n-1))
            omea = Rez(1,1);
        end
        if omea >= om_mx(gm)
            if om_mx(gm) > 1.005*G.ome(n-1)
                omea = om_mx(gm);
            else
                omea = 1.005*G.ome(n-1);
            end
        end
        if ~isempty(Rez)&&(omea == Rez(1,1))
            G1 = Rez(1,2);
        else
            if Tm > 0
                G1 = evalfr(Name,exp(1i*omea*Tm)); %feval(ProbA,p);
            else
                G1 = evalfr(Name,1i*omea);
            end
        end
        Kuv = 1;
        if (20*log10(abs(G1)) > -28) || (G.mdb(n-1) > -28)    %-28
             while 1
                 Mdb = [G.mdb(n-2), G.mdb(n-1), 20*log10(abs(G1))];
                 Ang = [G.ang(n-2), G.ang(n-1), angle(G1)];
                 [gam,sh] = crivizna_f(Mdb,Ang);              %crivizna_f
                 shd = dls + dls*((Mdb(3) + nidA)/nidA)^2;
                 if (gam < gam1)&&(sh < shd/4)
                     Kuv = 2;
                     break
                 end
                 if (gam < gam2)&&(sh < shd)
                     break
                 else    %gam >= gam2
                    Rez = [omea, G1; Rez];     %Засылка в Rez
                    omea = G.ome(n-1) + 0.6*(omea-G.ome(n-1));     %0.6!!!!!!!!!!!!!%%%Снижение omea
                    if Tm > 0
                        G1 = evalfr(Name, exp(1i*omea*Tm));%feval(ProbA,p);
                    else
                        G1 = evalfr(Name,1i*omea);
                    end
                 end    %gam >= gam2
             end        %while
        end    % if 20*log10(abs(G1))      %-28
        G.mdb(n) = 20*log10(abs(G1));   %Gout(n,3)
        G.ang(n) = angle(G1);   %Gout(n,4)
        ome = omea;
        G.ome(n) = ome;
        while (~isempty(Rez))&&(ome >= Rez(1,1))       %Чистка резерва
            Rez(1,:) = [];
        end
        if ome >= om_mx(gm)
            Nmx(gm) = n;
            gm = gm+1;
        end
        n = n + 1;            
    end    %n > 2
end     %while ome <= omk
disp('n')
disp(n-1)

            

Pzon = 0;  Vzon = []; XA = []; XF = []; 
Gout = [(1:n-1)',G.ome',G.mdb',G.ang'];

for n = 2:length(G.mdb)
                    %% Алгоритм нахожд. особых точек
    if (Gout(n,3)*Gout(n-1,3) <= 0)&&(abs(Gout(n,4)) < pi/2)   
       rab1 = abs(Gout(n-1,3)/(Gout(n,3)-Gout(n-1,3)));
       rab2 = Gout(n-1,1)+(Gout(n,1)-Gout(n-1,1))*rab1;  % N
       rab3 = Gout(n-1,2)+(Gout(n,2)-Gout(n-1,2))*rab1;  % om 
       rab4 = Gout(n-1,4)+(Gout(n,4)-Gout(n-1,4))*rab1;  % F
       %home
%disp('   И д ё т  р а с ч ё т  А Ф Ч Х.'   )
       XF = [XF; rab3 rab4*180/pi rab2];       % om F N
%        disp(['N=',int2str(Gout(n-1,1))]);
    end
    if (Gout(n,4)*Gout(n-1,4)<=0)&&(abs(Gout(n,4))<pi/2)...
                                     &&(abs(Gout(n-1,4))<pi/2)   
          rab1=abs(Gout(n-1,4)/(Gout(n,4)-Gout(n-1,4)));
          rab2=Gout(n-1,1)+(Gout(n,1)-Gout(n-1,1))*rab1;
          rab3=Gout(n-1,2)+(Gout(n,2)-Gout(n-1,2))*rab1;
          rab4=Gout(n-1,3)+(Gout(n,3)-Gout(n-1,3))*rab1;
        if rab4>=-nidA
         % home
            % disp('                 И д ё т  р а с ч ё т  А Ф Ч Х.')
            XA = [XA; rab3 rab4 rab2];           % om A N
%            disp(['N=',int2str(Gout(n-1,1))]);
        end
    end                    %% Алгоритм нахожд. особых точек 
                           %% Алгоритм контроля вх. в зону
    if abs(Gout(n,4)) <= pi/6
       rab1 = Rzon*sqrt(1 - (Gout(n,4)/(pi/6))^2);
       if (Gout(n,3) >= -rab1)&&(Gout(n,3) <= rab1)
          if Pzon == 0
             Vzon = [Vzon n-1];
             Pzon = 1;
          end
       else
          if Pzon == 1
             Vzon = [Vzon n];
             Pzon = 0;
          end
       end
    else
          if Pzon == 1
             Vzon = [Vzon n];
             Pzon = 0;
          end
    end                    
end%n

if ~isempty(Tunst)
    disp('Частоты (1/с) и постоянные времени (с) неустойчивых мод:')
    disp(Tunst)
end

              %%%%  	 ГРАФИКИ   	%%%%

Dat=fix(clock);
%%%%%%Обработка данных для построения графиков%%%
	% Опр. масшт. коэфф. при А 
rab1=-100; rab2=100;                    
for j=1:n                            
   if Gout(j,3)>rab1
      rab1=Gout(j,3);
   end
   if Gout(j,3)<rab2
      rab2=Gout(j,3);
   end
end
if abs(rab1)>abs(rab2)
   kga=abs(rab1);
else
   kga=abs(rab2);
end
if kga<1
    kga=1;
end

rab1=log10(200/kga)-floor(log10(200/kga));
kga=Vnorm(floor(rab1*10)+1)*10^floor(log10(200/kga));
if kga<4
   kga=4;
end         % Опр. масшт. коэфф. при А                         
%------------------------------------------------------------------------------------------%
switch nargin
     case 1, Ngraf=[1,3];
     case 2, Ngraf=gr1;
     case 3, Ngraf=[gr1,gr2];
     otherwise, Ngraf=1:3;
end
%-------------------------------График 1-----------------------------------------------------------%
 if any(Ngraf==1)                      
    figure%(1)
      Nf(1)=gcf;
      ln=log10(omn);
      lk=log10(omk);
      rab1=floor(ln);   X1=10^rab1; % Опр. масшт. по оси Х   
      rab2=ceil(lk);    X2=10^rab2; % Опр. масшт. по оси Х
     % Xmem(1)=X2*(1-0.25*(log10(X2)-log10(X1)));
      Xmem(1)=X2*(X1/X2)^0.2;
              % Основной график
      semilogx(Gout(:,2),kga*Gout(:,3),'b',Gout(:,2),Gout(:,4)*180/pi,'m+','Markersize',3);
    axis([X1 X2 -200 200]);		axis(axis);
     xlabel('Кругов. частота, 1/с'); 
     ylabel([num2str(kga) '*A, дб,  F, град']); grid;
      hold on
      if sum(size(XF))>0
         for j=1:length(XF(:,2))
            semilogx([XF(j,1) XF(j,1)],[200 XF(j,2)],'b-.')
            if j==1 
               text(XF(j,1)*(1+0.015*(rab2-rab1)),190,num2str(round(XF(j,2))));
            else
               if XF(j,1)/XF(j-1,1)-1>0.115*(rab2-rab1)  % 9mm
                  text(XF(j,1)*(1+0.015*(rab2-rab1)),190,num2str(round(XF(j,2)))); 
               end
            end
         end
      end
      if sum(size(XA))>0
         for j=1:length(XA(:,2))
            semilogx([XA(j,1) XA(j,1)],[-200 kga*XA(j,2)],'r--')
            if j==1 
               text(XA(j,1)*(1+0.015*(rab2-rab1)),-190,...
                        num2str(round(XA(j,2)*10)/10));
            else
               if XA(j,1)/XA(j-1,1)-1>0.17*(rab2-rab1)   % 14mm
                   text(XA(j,1)*(1+0.015*(rab2-rab1)),-190,...
                        num2str(round(XA(j,2)*10)/10));
               end
            end
         end
      end
    for j=1:2:(length(Vzon)-1)   % Прямоугольный контур 
         rab1=[Gout(Vzon(j),2) Gout(Vzon(j),2) Gout(Vzon(j+1),2) ...
                Gout(Vzon(j+1),2)];
         rab2=[130 150 150 130];
         semilogx(rab1,rab2,'r')
        semilogx(rab1,-rab2,'r')
    end
    if ~isempty(Tunst)%
       text(X1*(X2/X1)^0.85,-240,'Unstab ')
    end
    rab2=[int2str(Dat(1)),' ',int2str(Dat(2)),' ',int2str(Dat(3))];
    %text(X1*(1+0.2*(log10(X2)-log10(X1))),-240,[rab1,rab2])
    text(X1*(X2/X1)^0.05,-240,rab2)
    title(tit)
    hold off
    axis('normal');
end	%Ngraf==1
%-------------------------------График 2-----------------------------------------------------------%
 if any(Ngraf==2)
    figure%(2)	
    Nf(2)=gcf;
    rab1=log10(n)-floor(log10(n));    % Опр. масшт. по оси Х
    rab1=Vnorm(ceil(rab1*10)+1)*10^floor(log10(n)); X1=rab1;
    axis([0 rab1 -200 200]);	 axis(axis);
    Xmem(2)=0.8*X1;
    % Основной график
    plot(Gout(:,1),kga*Gout(:,3),'b',Gout(:,1),Gout(:,4)*180/pi,'m+',...
       Gout(:,1),100*log10(Gout(:,2)),'k-.','Markersize',3);
    axis([0 rab1 -200 200]);	 axis(axis);     
   xlabel('Номер точки N'); 
   ylabel([num2str(kga) '*A, дб,  F, град, 100*log10(om)']); grid;
    hold on
    if sum(size(XF))>0
     for j=1:length(XF(:,2))
        plot([XF(j,3) XF(j,3)],[200 XF(j,2)],'b-.')
        if (j==1)
           text(XF(j,3)+0.006*rab1,190,num2str(round(XF(j,2))));
        else
           if XF(j,3)-XF(j-1,3)>0.044*rab1      % 9mm
              text(XF(j,3)+0.006*rab1,190,num2str(round(XF(j,2))));  
           end
        end
     end
    end
    if sum(size(XA))>0
        for j=1:length(XA(:,2))
            plot([XA(j,3) XA(j,3)],[-200 kga*XA(j,2)],'r--')
            if (j==1)
               text(XA(j,3)+0.006*rab1,-190,num2str(round(XA(j,2)*10)/10));
            else
               if XA(j,3)-XA(j-1,3)>0.066*rab1      % 14mm
                  text(XA(j,3)+0.006*rab1,-190,num2str(round(XA(j,2)*10)/10));
               end
            end
        end
    end
    for j=1:2:(length(Vzon)-1)   % Прямоугольный контур 
         rab1=[Gout(Vzon(j),1) Gout(Vzon(j),1) Gout(Vzon(j+1),1) ...
                Gout(Vzon(j+1),1)];
         rab2=[130 150 150 130];
         plot(rab1,rab2,'r')
         plot(rab1,-rab2,'r')
    end
    if ~isempty(Tunst)%
        text(X1*8/10,-240,'Unstab')
    end
    rab2=[int2str(Dat(1)) ' ' int2str(Dat(2)) ' ' int2str(Dat(3))];
    text(X1/10,-240,rab2)
    title(tit)
    hold off
    axis('normal');
end	%Ngraf==2
%-------------------------------График 3-----------------------------------------------------------%
if any(Ngraf==3)
    figure%(3)	
%     Vaxis=[-90 40 -50 50];
    Vaxis=[-100 30 -50 50];
    Nf(3)=gcf;
    for j=1:n
        if Gout(j,3)>-nidA		% Gout(j,3) - 20*log10(A)
           rab1=Gout(j,3)+nidA;
        else
           rab1=0;
        end
        Rgod(j)=rab1;
        XsmY(j,1)=rab1*cos(Gout(j,4))-nidA;     %Координаты X,Y
        XsmY(j,2)=rab1*sin(Gout(j,4));
    end
	

      rab1=pi/12;         % Рисование сетки
      for j=1:25
         Ring(j,1:2)=[cos(j*rab1)  sin(j*rab1)];
      end
      rab1=nidA*ones(size(Ring(:,1)));
      rab2=(rab1+Rzon*Ring(:,1)).*cos(pi/6*Ring(:,2));
      rab3=(rab1+Rzon*Ring(:,1)).*sin(pi/6*Ring(:,2));
      axis(Vaxis);
	  axis(axis);
      hold on
  
      plot([Vaxis(1),Vaxis(2),Vaxis(2)],[Vaxis(4),Vaxis(4),Vaxis(3)],'k')   % Рамка
      patch(rab2-nidA,rab3,[0.9,0.9,0.9],'EdgeColor','none')	% Контур окрестности крит. точки
      patch((rab2-nidA)*2/3,rab3*2/3,'w','EdgeColor','none')	% Контур окрестности крит. точки

      rab1=[-90,-80,-70,-50:10:-10,10,20]; rab2=zeros(size(rab1));		% Метки 10 дб +
      plot(rab1, rab2, '+k','Markersize',4)  %горизонтальные
	  rab1=-20:10:20; rab2=-nidA*ones(size(rab1));	
      plot(rab2, rab1, '+k','Markersize',4)  %вертикальные
      
      plot((Ring(1:23,1)-1)*nidA,Ring(1:23,2)*nidA,'xk','Markersize',5); % Окружн. един. радиуса
      plot(0, 0, '+','Markersize',12);        plot(0, 0, 'o')   %Метка в крит. точке

        %%%%%%%%%%%%%%%%%%%% Основной график %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Col=['r  ';'b-o';'k--';'m-*';'c--';'g  ';'m-s'];
      Colar=['r';'b';'k';'m';'c';'y';'g'];	Ncl=7;
		s=0;   Vr=XsmY(1,:)';
      scol=rem(s,Ncl)+1;
for j=2:n
    rab1=~(((Rgod(j-1)<10)&&(Rgod(j)>=10))||(j==n));
    if rab1	% Непересечение -20дб вверх
       Vr=[Vr,XsmY(j,:)'];
    else
       Vr=[Vr,XsmY(j,:)'];
       scol=rem(s,Ncl)+1;
     %plot(Vr(1,:),Vr(2,:),'r'); 
     plot(Vr(1,:),Vr(2,:),Col(scol,:),'Markersize',3); % Рисование линии текущего цвета
       s=s+1;
       Vr=XsmY(j,:)';
    end

    if (Gout(j,2)-1)*(Gout(j-1,2)-1)<=0	%Переход через om=1
      rab=((XsmY(j,1)-XsmY(j-1,1))^2+(XsmY(j,2)-XsmY(j-1,2))^2)^0.5;
       Ear=((XsmY(j,1)-XsmY(j-1,1))+i*(XsmY(j,2)-XsmY(j-1,2)))/rab;
            Rotarr=(arrow*Ear).';
       patch(XsmY(j,1)+real(Rotarr),XsmY(j,2)+imag(Rotarr),'r');	% Стрелка om=1
    end
end
%--------------------------Частоты на выносках----------------------------
Mt1=-64: 12:-16; %Mt2=43*ones(1,5);	% Массив координат выносок значений частот
rab=Rgod(Nmx);
for s=(length(rab)-1):-1:1
    if rab(s)==0
        Nmx(s)=[];
    end
end
if length(Nmx)>10
    Nmx=Nmx(1:10);
end
Ctf=[XsmY(Nmx,:),Gout(Nmx,2)];
[rab,I]=sort(Ctf(:,2));
Ctf=Ctf(I,:);
Ctfo=Ctf(1:fix(length(Nmx)/2),:);
Ctfp=Ctf(fix(length(Nmx)/2)+1:length(Nmx),:);
[rab,I]=sort(Ctfo(:,1));
Ctfo=Ctfo(I,:);
[rab,I]=sort(Ctfp(:,1));
Ctfp=Ctfp(I,:);
nnao=fix((5-length(Ctfo(:,1)))/2);
nnap=fix((5-length(Ctfp(:,1)))/2);
for s=1:length(Ctfp(:,1))
    plot(Ctfp(s,1),Ctfp(s,2),'.')
    plot([Ctfp(s,1),Mt1(s+nnap)],[Ctfp(s,2),43],':g')
    plot([Mt1(s+nnap),Mt1(s+nnap)+6],[43,43],'linew',2)
    text(Mt1(s+nnap),46,num2str(round(Ctfp(s,3)*10)/10));
end
for s=1:length(Ctfo(:,1))
    plot(Ctfo(s,1),Ctfo(s,2),'.')
    plot([Ctfo(s,1),Mt1(s+nnap)],[Ctfo(s,2),-43],':g')
    plot([Mt1(s+nnao),Mt1(s+nnao)+6],[-43,-43],'linew',2)
    text(Mt1(s+nnap),-46,num2str(round(Ctfo(s,3)*10)/10));
end

 text(16,-41,[int2str(PZ),'    ',int2str(Nppo)]);
 text(15,-46,[' T=',num2str(Tm)]);
      if ~isempty(Tunst)%
        text(16,-36,'Unstab')
      end
      if sum(size(XA))>0
         XAg=[round(XA(:,1)*10)/10 round(XA(:,2)*10)/10];
         text(3,43,'om         A')
         for j=1:length(XA(:,2))
            text(3,42-j*4,num2str(XAg(j,1)))
            text(18,42-j*4,num2str(XAg(j,2)))
         end
      end
      if sum(size(XF))>0
         XFg=[round(XF(:,1)*10)/10 round(XF(:,2))];
         text(-94,43,'om')
         text(-78,43,'F')
         for j=1:length(XF(:,2))
            text(-94,42-j*4,num2str(XFg(j,1)))
            text(-79,42-j*4,num2str(XFg(j,2)))
         end
      end
      text(-94,-46,[int2str(Dat(1)) ' ' int2str(Dat(2)) ' ' int2str(Dat(3))])
      title(tit);
      hold off
      set(gca,'XTick',[],'YTick',[])
      axis('normal');
 end	%Ngraf==3


function [gam,sh]=crivizna_f(Mdb,Ang)
Mdb=Mdb+30;
for s=1:3
    if Mdb(s)<0
        Mdb(s)=0;
    end
end
a=sqrt(Mdb(1)^2+Mdb(2)^2-2*Mdb(1)*Mdb(2)*cos(Ang(1)-Ang(2)));
b=sqrt(Mdb(3)^2+Mdb(2)^2-2*Mdb(3)*Mdb(2)*cos(Ang(3)-Ang(2)));
c=sqrt(Mdb(1)^2+Mdb(3)^2-2*Mdb(1)*Mdb(3)*cos(Ang(1)-Ang(3)));
p=sum([a b c])/2;
r=sqrt((p-a)*(p-b)*(p-c)/p);
if (p-b)>eps%(a>eps)&&(b>eps)&&(c>eps)
    gam=2*r/(p-b);
else
    gam=0.02;
end
sh=b;

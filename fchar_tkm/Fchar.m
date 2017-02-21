%FCHAR5 - ѕрограмма расчЄта частотных характеристик систем (аналоговых 
%и дискретных с одной частотой прерывания)
clc
disp('       ѕрограмма расчЄта частотных характеристик дискретной системы ')
disp('                   с одной частотой прерывания FCHAR')
clear all
clear global
global  INQ TIN
%i=(-1)^0.5;
dlr=6;         % Ќом. расстояние между точками полярн. графика, дб
obA=20;        % ¬ерхний уровень амплитуд полярн. графика, дб
nidA=30;       % ”ровень нуля амплитуд полярн. графика, дб
Rzon=6;        % –адиус окрестности крит. точки
Vnorm=[1 1.25 1.6 2 2.5 3 4 5 6 8 10];  % –яд чисел для границ шкалы графиков
ProbD=' ';
rab=3^0.5/2;
arrow=2*[-rab+i/2;0;-rab-i/2;-rab+i/2];
%------------------------------------------------------------------------------------------%
load fcharust
disp(' ')
disp([blanks(10),'ѕредыдущие установки fchar: omn=',num2str(omn)])
disp(['                                      omk=',num2str(omk)])
disp(['                                       Nd=',num2str(Nd)])
rab1=input('ћеняем? (y/<Enter>): ','s');
disp(' ')
if strcmp(rab1,'y')
	rab1=0;
	while rab1==0
        	omn=input('¬ведите начальную частоту omn=');
        	omk=input('¬ведите  конечную частоту omk=');
        	if omk>omn 
                	rab1=1;
        	else
           	disp(' ')
           	disp('ќшибка ввода. ¬ведите omk>omn');
        	end
	end
	Nd=input('¬ведите минимальное число шагов на декаду Nd='); 
	disp(' ')
end
save fcharust omn omk Nd -append
disp([blanks(10),'ѕредыдущие установки fchar: '])
disp(['- имя m-файла  аналоговой структуры: ',ProbA])
disp(' ')
disp('¬ведите новое имя m-файла описания аналоговой структуры ');
ProbAin=input('(для сохранения прежнего - <Enter>): ','s');
disp(' ')
if length(ProbAin)>0 ProbA=ProbAin; end
tit=ProbA;
INQ=1;         % «апрос диалога из рабочего файла
C1=feval(ProbA,i*omn);
disp([blanks(10),'ѕредыдущие установки fchar: TIN=',num2str(TIN)])
 	disp(' ')
	disp('¬ведите новую величину периода прерывания TIN') 
	disp('            (для аналоговой системы TIN=0);')
	Tinin=input('(для сохранения прежнего - <Enter>): '); 
if length(Tinin)>0 TIN=Tinin; end
if TIN>0
   disp(' ')
	disp([blanks(10),'ѕредыдущие установки fchar: '])
	disp(['- имя m-файла дискретной структуры: ',ProbD])
	disp(' ')
   disp('¬ведите новое имя m-файла описания дискретной структуры ');
   ProbDin=input('(для сохранения прежнего - <Enter>): ','s');
   disp(' ')
   if length(ProbDin)>0 ProbD=ProbDin; end
	tit=[tit '  ' ProbD '  T=' num2str(TIN)];
        C1=feval(ProbD,C1,i*omn,TIN);      % «апрос диалога из рабочего файла
        if omk>pi/TIN
           omk=pi/TIN;
           disp('«начение omk ограничено половиной частоты прерывания')
           disp(' ')
        end
end
INQ=0;         % «апрос диалога из рабочего файла отключен
disp(' ')
disp([blanks(10),'ѕредыдущие установки fchar: DTEXT=',DTEXT])
disp(' ')
rab1=input('ћеняем DTEXT? (y/<Enter>): ','s');
disp(' ')
if strcmp(rab1,'y')
   disp('¬ведите новый дополнительный текст в заголовок графика DTEXT')
   DTEXT=input('(для отсутствия - <Enter>): ','s');
   disp(' ')
end
if length(DTEXT)>0    tit=[tit '  ' DTEXT]; end
save  fcharust  ProbA TIN ProbD DTEXT  -append
disp(' ')
%clc
disp('                 » д Є т  р а с ч Є т  ј ‘ „ ’.')
disp(' ')
ln=log10(omn); lt=ln; lk=log10(omk); i=sqrt(-1); nfc=1; Ndt=Nd;
XA=[]; XF=[]; Gout=[]; Pzon=0; Vzon=[]; SumAng=0; Prg=1;

while lt<=lk
  Gout(nfc,1)=nfc;                % номер шага по частоте
  Gout(nfc,2)=10^lt;            % omega
  p=i*Gout(nfc,2);
  C1=feval(ProbA,p);
  if TIN>0
       si=-1;
       Cin=C1/TIN;
       omsum=Gout(nfc,2)+si*2*pi/TIN;
       while abs(omsum)<8/TIN         % Ўирина полосы частот
          Cin=Cin+feval(ProbA,i*omsum)/TIN;
          si=si-1;
          omsum=Gout(nfc,2)+si*2*pi/TIN;
       end
       si=1;
       omsum=Gout(nfc,2)+si*2*pi/TIN;
       while abs(omsum)<8/TIN         % Ўирина полосы частот
          Psum=i*omsum;
          Cin=Cin+feval(ProbA,i*omsum)/TIN;
          si=si+1;
          omsum=Gout(nfc,2)+si*2*pi/TIN;
       end
       z=exp(p*TIN);
       C1=feval(ProbD,Cin,z,TIN);
  end
  if length(C1)>1
     error('ќшибка описания системы. ¬ыход должен быть скаляром.');
  end
  Gout(nfc,3)=20*log10(abs(C1));
  Gout(nfc,4)=angle(C1);
 if abs(C1-1)>1e-10
     Gout(nfc,5)=angle(C1-1);
  else
     Gout(nfc,5)=0;
  end
% Gout: [N,om,A,F(rad),Fsumangl(rad)]
  if nfc>1                          
%%%%%%%%%%%%%%%%%%%%%% јлгоритм изменения величины шага %%%%%%%%%%%%%%%%%%%%%%%%%%%
     Rad=[Gout(nfc,3) Gout(nfc-1,3)];
     nidAst=nidA+10;
     for j=1:2
       if Rad(j)>-nidAst
         Rad(j)=Rad(j)+nidAst;
       else
         Rad(j)=0;
       end
     end
     if Rad(1)>nidAst			%% dlrc - требуемая текущая dl
	     dlrc=dlr+2*dlr*Gout(nfc,3)/nidA;
     else
	     dlrc=dlr/(1+2*(nidAst-Rad(1))/nidAst);	% ¬ нуле dlrc уменьшен втрое
     end	
     dlrt=sqrt(Rad*Rad'-2*Rad(1)*Rad(2)*cos(Gout(nfc,4)-Gout(nfc-1,4)));
     if (dlrt<1.5*dlrc)|(Ndt>2000)  
         if (dlrt<1.5*dlrc) %Ўаг годится
            Ndt=round(Ndt*dlrt/dlrc);	%% dlrc       
            if Ndt<Nd
               Ndt=Nd;
            end
         else  %Ndt>2000
            disp('¬нимание! ќсобая точка на мнимой оси.')
            disp(['Omega=',num2str(Gout(nfc,2))])            
         end
        %-------------------% јлгоритм нахожд. точек A=0, F=0 %--------------------%
                 % Gout: [N,om,A,F(rad),Fsumangl(rad)]
         if Gout(nfc,3)>=-(nidA+10)
            if (Gout(nfc,3)*Gout(nfc-1,3)<=0)&(abs(Gout(nfc,4))<pi/2)   
               rab1=abs(Gout(nfc-1,3)/(Gout(nfc,3)-Gout(nfc-1,3)));
               rab2=Gout(nfc-1,1)+(Gout(nfc,1)-Gout(nfc-1,1))*rab1;  % N
               rab3=Gout(nfc-1,2)+(Gout(nfc,2)-Gout(nfc-1,2))*rab1;  % om 
               rab4=Gout(nfc-1,4)+(Gout(nfc,4)-Gout(nfc-1,4))*rab1;  % F
               %home
               XF=[XF; rab3 rab4*180/pi rab2];       % om F N
               disp(['N=',int2str(Gout(nfc-1,1))])
            end
            if (Gout(nfc,4)*Gout(nfc-1,4)<=0)&(abs(Gout(nfc,4))<pi/2)...
                                         &(abs(Gout(nfc-1,4))<pi/2)   
                rab1=abs(Gout(nfc-1,4)/(Gout(nfc,4)-Gout(nfc-1,4)));
                rab2=Gout(nfc-1,1)+(Gout(nfc,1)-Gout(nfc-1,1))*rab1;
                rab3=Gout(nfc-1,2)+(Gout(nfc,2)-Gout(nfc-1,2))*rab1;
                rab4=Gout(nfc-1,3)+(Gout(nfc,3)-Gout(nfc-1,3))*rab1;	% ј
                if rab4>=-30
                  %home
                  XA=[XA; rab3 rab4 rab2];           % om A N
                  disp(['N=',int2str(Gout(nfc-1,1))])
               end
            end   %% јлгоритм нахожд. точек A=0, F=0 
            %---------------% јлгоритм контроля вх. в зону %-----------------------%
            if abs(Gout(nfc,4))<=pi/6
               rab1=Rzon*sqrt(1-(Gout(nfc,4)/(pi/6))^2);
               if (Gout(nfc,3)>=-rab1)&(Gout(nfc,3)<=rab1)
                  if Pzon==0
                     Vzon=[Vzon nfc-1];
                     Pzon=1;
                  end
               else
                  if Pzon==1
                     Vzon=[Vzon nfc];
                     Pzon=0;
                  end
               end
            else
                  if Pzon==1
                     Vzon=[Vzon nfc];
                     Pzon=0;
                  end
            end  %% јлгоритм контроля вх. в зону            
         end %if Gout(nfc,3)>=-(Nida+10)
         %---------------------------------------------------------------------------%
        rab1=Gout(nfc,5)-Gout(nfc-1,5);
        if abs(rab1)>pi
           rab1=rab1-2*pi*sign(rab1);   
        end
        SumAng=SumAng+rab1;
        lt=lt+1/Ndt;
        nfc=nfc+1;
     else   % (if dlrt>=1.5*dlrc)&(Ndt<=2000), шаг негоден, возврат
        lt=lt-1/Ndt;
        Ndt=round(Ndt*dlrt/dlrc);
        lt=lt+1/Ndt;
     end     %if dlrt<1.5*dlrc   % јлгоритм изменения величины шага
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  else                                      % nfc=1
     lt=lt+1/Ndt;
     nfc=nfc+1;
  end
end
%------------------------------------------------------------------------------------------%
nfc=nfc-1;
if Pzon==1
   Vzon=[Vzon nfc];
   Pzon=0;
end
SumAng=round(SumAng*10)/10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%  	 √–ј‘» »   	%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dat=fix(clock);
disp('         ќбработка данных для построения графиков.')
	% ќпр. масшт. коэфф. при ј 
      rab1=-100; rab2=100;                    
      for j=1:nfc                            
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
      end         % ќпр. масшт. коэфф. при ј                         

   %------------------------------------------------------------------------------------------%         
 figure%(1)
 axis('normal');
      Nf(1)=gcf;
      rab1=floor(ln);   X1=10^rab1; % ќпр. масшт. по оси ’   
      rab2=ceil(lk);    X2=10^rab2; % ќпр. масшт. по оси ’
     % Xmem(1)=X2*(1-0.25*(log10(X2)-log10(X1)));
      Xmem(1)=X2*(X1/X2)^0.2;
              % ќсновной график
      semilogx(Gout(:,2),kga*Gout(:,3),'b',Gout(:,2),Gout(:,4)*180/pi,'m+','Markersize',3);
  	axis([X1 X2 -200 200]);		axis(axis);
     xlabel(' ругов. частота, 1/с'); 
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
      for j=1:2:(length(Vzon)-1)   % ѕрямоугольный контур 
         rab1=[Gout(Vzon(j),2) Gout(Vzon(j),2) Gout(Vzon(j+1),2) ...
                Gout(Vzon(j+1),2)];
         rab2=[130 150 150 130];
         semilogx(rab1,rab2,'r')
      semilogx(rab1,-rab2,'r')
      end
     	rab1=[num2str(SumAng),'  ',int2str(Dat(1)),' '];
	rab2=[int2str(Dat(2)),' ',int2str(Dat(3))];
      %text(X1*(1+0.2*(log10(X2)-log10(X1))),-240,[rab1,rab2])
      text(X1*(X2/X1)^0.05,-240,[rab1,rab2])
      title(tit)
	hold off
   %------------------------------------------------------------------------------------------%
figure%(2)	
axis('normal');
      Nf(2)=gcf;
      rab1=log10(nfc)-floor(log10(nfc));    % ќпр. масшт. по оси ’
      rab1=Vnorm(ceil(rab1*10)+1)*10^floor(log10(nfc)); X1=rab1;
      axis([0 rab1 -200 200]);	 axis(axis);
      Xmem(2)=0.8*X1;
	 % ќсновной график
      plot(Gout(:,1),kga*Gout(:,3),'b',Gout(:,1),Gout(:,4)*180/pi,'m+',...
           Gout(:,1),100*log10(Gout(:,2)),'k-.','Markersize',3);
	axis([0 rab1 -200 200]);	 axis(axis);     
   xlabel('Ќомер точки N'); 
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
      for j=1:2:(length(Vzon)-1)   % ѕрямоугольный контур 
         rab1=[Gout(Vzon(j),1) Gout(Vzon(j),1) Gout(Vzon(j+1),1) ...
                Gout(Vzon(j+1),1)];
         rab2=[130 150 150 130];
         plot(rab1,rab2,'r')
         plot(rab1,-rab2,'r')
      end
      	rab1=[num2str(SumAng),'  ',int2str(Dat(1)),' '];
	rab2=[int2str(Dat(2)),' ',int2str(Dat(3))];
      text(X1/10,-240,[rab1,rab2])
      title(tit)
	hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure%(3)	
axis('normal');
Vaxis=[-90 40 -50 50];
       Nf(3)=gcf;
% Gout: [N,om,A,F(rad),Fsumangl(rad)]
for j=1:nfc
   if Gout(j,3)+nidA>0		% Gout(j,3) -   20*log10(A)
           rab1=Gout(j,3)+nidA;
   else
           rab1=0;
   end
	Rgod(j)=rab1;			% –адиус
   XsmY(j,1)=rab1*cos(Gout(j,4))-nidA;  % -ты графика (декарт. к-ты 
   XsmY(j,2)=rab1*sin(Gout(j,4));       %относительно крит. точки)
   %------------------------------------------------------------------------------------------%
                     %»сключение перечеркиваний
   if  XsmY(j,1)<Vaxis(1)
      XsmY(j,1)=Vaxis(1);
   elseif XsmY(j,1)>Vaxis(2)
      XsmY(j,1)=Vaxis(2);
   end
   if  XsmY(j,2)<Vaxis(3)
      XsmY(j,2)=Vaxis(3);
   elseif XsmY(j,2)>Vaxis(4)
      XsmY(j,2)=Vaxis(4);
   end
end
%------------------------------------------------------------------------------------------%
rab1=pi/12;         % –исование сетки
for j=1:25
         Ring(j,1:2)=[cos(j*rab1)  sin(j*rab1)];
end
      rab1=nidA*ones(size(Ring(:,1)));
      rab2=(rab1+Rzon*Ring(:,1)).*cos(pi/6*Ring(:,2));
      rab3=(rab1+Rzon*Ring(:,1)).*sin(pi/6*Ring(:,2));
      axis(Vaxis);
	  axis(axis);
      hold on
      
      plot([Vaxis(1),Vaxis(2),Vaxis(2)],[Vaxis(4),Vaxis(4),Vaxis(3)],'k')   % –амка
      patch(rab2-nidA,rab3,[0.9,0.9,0.9],'EdgeColor','none')	%  онтур окрестности крит. точки
      patch((rab2-nidA)*2/3,rab3*2/3,'w','EdgeColor','none')	%  онтур окрестности крит. точки
      
      rab1=[-80,-70,-50:10:-10,10,20,30]; rab2=zeros(size(rab1));		% ћетки 10 дб +
      plot(rab1, rab2, '+k','Markersize',4)  %горизонтальные
	  rab1=-20:10:20; rab2=-nidA*ones(size(rab1));	
      plot(rab2, rab1, '+k','Markersize',4)  %вертикальные
      
      plot((Ring(1:23,1)-1)*nidA,Ring(1:23,2)*nidA,'xk','Markersize',5); % ќкружн. един. радиуса
      plot(0, 0, '+','Markersize',12);        plot(0, 0, 'o')   %ћетка в крит. точке
      
      %%%%%%%%%%%%%%%%%%%% ќсновной график %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      Mt1=-54: 12:-6; Mt2=43*ones(1,5);
      Mt=[Mt1,Mt1; Mt2,-Mt2];	% ћассив координат выносок значений частот
      Col=['r  ';'k: ';'m* ';'c--';'b: ';'y  ';'g--'];
      Colar=['r';'k';'m';'c';'b';'y';'g'];	Ncl=7;
		s=0;   Vr=XsmY(1,:)';
for j=2:nfc
   %rab1=~(((Rgod(j-1)<10)&(Rgod(j)>=10))|(j==nfc));
   %if rab1	% Ќепересечение -20дб вверх
   if j>2  
      rab1=(Rgod(j-1)-Rgod(j-2))*(Rgod(j)-Rgod(j-1));
   else   
      rab1=1;
   end
   if ((rab1<0)&(Rgod(j-1)-Rgod(j-2)>0)&(Rgod(j)>10))|(j==nfc)	%ћаксимум радиуса
	     Vr=[Vr,XsmY(j,:)'];
	     scol=rem(s,Ncl)+1;
         %plot(Vr(1,:),Vr(2,:),'r'); 
         plot(Vr(1,:),Vr(2,:),Col(scol,:),'Markersize',3); % –исование линии текущего цвета
	     s=s+1;
	     Vr=XsmY(j,:)';
        %------------------------------------------------------------------------------------% 
     if Rgod(j)>10
         plot(XsmY(j,1),XsmY(j,2),'.')
                    % —трелка:
				% Gout: [N,om,A,F(rad),Fsumangl(rad)]
         rab=abs(Gout(j,4)-Gout(j-1,4))<pi;
         rab1=-(~rab-rab)*pi/2;        
         if (Gout(j,4)<Gout(j-1,4))
            Ear=cos(Gout(j,4)-rab1)+i*sin(Gout(j,4)-rab1);
         else
            Ear=cos(Gout(j,4)+rab1)+i*sin(Gout(j,4)+rab1);
         end
         Rotarr=(arrow*Ear).';
         patch(XsmY(j,1)+real(Rotarr),XsmY(j,2)+imag(Rotarr),Colar(scol,:));	% —трелка	   
 		   if ~isempty(Mt)
            [L,k]=min((XsmY(j,1)-Mt(1,:)).^2+(XsmY(j,2)-Mt(2,:)).^2);
            plot([XsmY(j,1),Mt(1,k)],[XsmY(j,2),Mt(2,k)],':m') %Ћинии выносок
            plot(Mt(1,k),Mt(2,k),'.')
            text(Mt(1,k),Mt(2,k)+3*sign(Mt(2,k)),num2str(round(Gout(j,2)*10)/10));
            rab=length(Mt(1,:));
            if rab==1
               Mt=[];
            else	%¬ырезание использованной выноски
               switch k
               case rab,Mt=Mt(:,1:k-1);
               case 1,Mt=Mt(:,2:rab);
               otherwise,Mt=[Mt(:,1:k-1),Mt(:,k+1:rab)];
               end
            end
         end
      end   %if Rgod(j)>10   
   else	%ћонотонность
	     Vr=[Vr,XsmY(j,:)'];
   end	%%ћаксимум радиуса
    %------------------------------------------------------------------------------------%     
   if (Gout(j,2)-1)*(Gout(j-1,2)-1)<=0	%ѕереход через om=1
      ((XsmY(j,1)-XsmY(j-1,1))^2+(XsmY(j,2)-XsmY(j-1,2))^2)^0.5;
       Ear=((XsmY(j,1)-XsmY(j-1,1))+i*(XsmY(j,2)-XsmY(j-1,2)))/ans;
            Rotarr=(arrow*Ear).';
       patch(XsmY(j,1)+real(Rotarr),XsmY(j,2)+imag(Rotarr),Colar(1,:));	% —трелка om=1
   end
end  %for j=2:nfc

 %------------------------------------------------------------------------------------%     

      text(-88,-46,['”ѕ = ',num2str(SumAng),' рад'])
      if sum(size(XA))>0
         XAg=[round(XA(:,1)*10)/10 round(XA(:,2)*10)/10];
         text(13,43,'om')
         text(29,43,'A')
         for j=1:length(XA(:,2))
            text(13,42-j*4,num2str(XAg(j,1)))
            text(28,42-j*4,num2str(XAg(j,2)))
         end
      end
      if sum(size(XF))>0
         XFg=[round(XF(:,1)*10)/10 round(XF(:,2))];
         text(-84,43,'om')
         text(-68,43,'F')
         for j=1:length(XF(:,2))
            text(-84,42-j*4,num2str(XFg(j,1)))
            text(-69,42-j*4,num2str(XFg(j,2)))
         end
      end
      text(20,-46,[int2str(Dat(1)) ' ' int2str(Dat(2)) ' ' int2str(Dat(3))])
      title(tit);
      hold off
      set(gca,'XTick',[],'YTick',[])    %»сключение оцифровки осей
disp(' ')
disp('ƒля сохранения графиков используйте функцию figsav (hgsave для версий MATLAB R11, R12).')      
disp(' ')
disp('ѕоздравляю с расчЄтом по программе FCHAR !')
      
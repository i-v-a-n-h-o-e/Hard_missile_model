function bode_os(Name,gr1,gr2,gr3)
%������ ������ � ������� �������� ������������ bode_os(Name,gr1,gr2,gr3) 
dlr=6;         % ���. ���������� ����� ������� ������. �������, ��
obA=20;        % ������� ������� �������� ������. �������, ��
nidA=30;       % ������� ���� �������� ������. �������, ��
Rzon=6;        % ������ ����������� ����. �����
Vnorm=[1 1.25 1.6 2 2.5 3 4 5 6 8 10];  % ��� ����� ��� ������ ����� ��������
Nd=20;         % ���. ����� ����� �� ������ �������
rab=3^0.5/2;
arrow=2*[-rab+i/2;0;-rab-i/2;-rab+i/2]; %������ ������� �������
%--------------------------------------------------------------------------%
Ts=Name.ts;    
if isempty(Name.notes)
   Nnotes.Notes=[];
else
   Nnotes=cell2struct(Name.notes,'Notes',1);
end
tit=[Nnotes.Notes,'  Ts=',num2str(Ts)];
%--------------------------------------------------------------------------%
Pol=pole(Name);
Nsys=length(Pol);
disp(['������� ������� ', Nnotes.Notes,' ����� ',num2str(Nsys)])
if Ts>0
   Npp=sum(abs(Pol)>1)
   omk=10^(floor(log10(pi/Ts)*10)/10);
   Ang=sort(abs(angle(Pol)));
   rab=0;
   s=1;
      while rab<=eps
       rab=Ang(s);
       s=s+1;
   end
   omn=1.2*rab/Ts;
else
   Npp=sum(real(Pol)>0)
   Im=sort(abs(imag(Pol)));
   omk=max(Im);
   rab=0;
   s=1;
   while rab<=eps
       rab=Im(s);
       s=s+1;
   end
   omn=rab;
end
if omk<100*omn
    rab=sqrt(omn*omk);
    omk=20*rab;
    omn=rab/10*2;
end
% omn
% omk
%--------------------------------------------------------------------------%
ln=log10(omn); lt=ln; lk=log10(omk);  nfc=1; Ndt=Nd;    %i=sqrt(-1);
XA=[]; XF=[]; Gout=[]; Pzon=0; Vzon=[]; SumAng=0; Prg=1;

while lt<=lk
  Gout(nfc,1)=nfc;                % ����� ���� �� �������
  Gout(nfc,2)=10^lt;            % omega
  p=i*Gout(nfc,2);
  if Ts>0
     C1=evalfr(Name,exp(p*Ts));%feval(ProbA,p);
  else
     C1=evalfr(Name,p);
  end
  if length(C1)>1
     error('������ �������� �������. ����� ������ ���� ��������.');
  end
  Gout(nfc,3)=20*log10(abs(C1));
  Gout(nfc,4)=angle(C1);
  if abs(C1-1)>1e-10
     Gout(nfc,5)=angle(C1-1);
  else
     Gout(nfc,5)=0;
  end

  if nfc>1                          
                  % �������� ��������� �������� ����
     Rad=[Gout(nfc,3) Gout(nfc-1,3)];
     nidAst=nidA+10;
     for j=1:2
       if Rad(j)>-nidAst
         Rad(j)=Rad(j)+nidAst;
       else
         Rad(j)=0;
       end
     end
     if Rad(1)>nidAst			%% ������ dlrc
	dlrc=dlr+2*dlr*Gout(nfc,3)/nidA;
     else
	dlrc=dlr/(1+2*(nidAst-Rad(1))/nidAst);	% � ���� dlrc �������� �����
     end	
     dlrt=sqrt(Rad*Rad'-2*Rad(1)*Rad(2)*cos(Gout(nfc,4)-Gout(nfc-1,4)));
     if dlrt<1.5*dlrc
 
       Ndt=round(Ndt*dlrt/dlrc);	%% dlrc
 
        if Ndt<Nd
          Ndt=Nd;
        end
                        %% �������� ������. ������ �����
        if (Gout(nfc,3)*Gout(nfc-1,3)<=0)&(abs(Gout(nfc,4))<pi/2)   
           rab1=abs(Gout(nfc-1,3)/(Gout(nfc,3)-Gout(nfc-1,3)));
           rab2=Gout(nfc-1,1)+(Gout(nfc,1)-Gout(nfc-1,1))*rab1;  % N
           rab3=Gout(nfc-1,2)+(Gout(nfc,2)-Gout(nfc-1,2))*rab1;  % om 
           rab4=Gout(nfc-1,4)+(Gout(nfc,4)-Gout(nfc-1,4))*rab1;  % F
           %home
  %disp('                 � � � �  � � � � � �  � � � �.')
           XF=[XF; rab3 rab4*180/pi rab2];       % om F N
%            disp(['N=',int2str(Gout(nfc-1,1))]);
        end
        if (Gout(nfc,4)*Gout(nfc-1,4)<=0)&(abs(Gout(nfc,4))<pi/2)...
                                         &(abs(Gout(nfc-1,4))<pi/2)   
              rab1=abs(Gout(nfc-1,4)/(Gout(nfc,4)-Gout(nfc-1,4)));
              rab2=Gout(nfc-1,1)+(Gout(nfc,1)-Gout(nfc-1,1))*rab1;
              rab3=Gout(nfc-1,2)+(Gout(nfc,2)-Gout(nfc-1,2))*rab1;
              rab4=Gout(nfc-1,3)+(Gout(nfc,3)-Gout(nfc-1,3))*rab1;
              if rab4>=-30
         %    home
  %disp('                 � � � �  � � � � � �  � � � �.')
            XA=[XA; rab3 rab4 rab2];           % om A N
%            disp(['N=',int2str(Gout(nfc-1,1))]);
            end
        end                    %% �������� ������. ������ ����� 
                               %% �������� �������� ��. � ����
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
        end                    %% �������� �������� ��. � ����
        rab1=Gout(nfc,5)-Gout(nfc-1,5);
        if abs(rab1)>pi
           rab1=rab1-2*pi*sign(rab1);   
        end
        SumAng=SumAng+rab1;
        lt=lt+1/Ndt;
        nfc=nfc+1;
     else
        %disp(['�������     ' num2str(nfc) '    ' num2str(Gout(nfc,2)) ...
        %        '     ' int2str(Ndt)]);
        lt=lt-1/Ndt;
        Ndt=round(Ndt*dlrt/dlrc);	%% dlrc
        lt=lt+1/Ndt;
     end                % �������� ��������� �������� ����
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
SumAng=round(SumAng/pi*10)/10;
disp(['SumAng=',num2str(SumAng),' pi'])
%------------------------------------------------------------------------------------------%
              %%%%  	 �������   	%%%%

Dat=fix(clock);
%disp('         ��������� ������ ��� ���������� ��������.')
	% ���. �����. �����. ��� � 
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
end         % ���. �����. �����. ��� �                         
%------------------------------------------------------------------------------------------%
Numgraf=nargin-1;
if Numgraf==0
   gr1=3;
   Numgraf=1;
end
for qgr=1:Numgraf
   switch qgr
     case 1, Ngraf=gr1;
     case 2, Ngraf=gr2;
     case 3, Ngraf=gr3;
     end
%------------------------------------------------------------------------------------------%
     if Ngraf==1                      
   figure%(1)
      Nf(1)=gcf;
      rab1=floor(ln);   X1=10^rab1; % ���. �����. �� ��� �   
      rab2=ceil(lk);    X2=10^rab2; % ���. �����. �� ��� �
     % Xmem(1)=X2*(1-0.25*(log10(X2)-log10(X1)));
      Xmem(1)=X2*(X1/X2)^0.2;
              % �������� ������
      semilogx(Gout(:,2),kga*Gout(:,3),'b',Gout(:,2),Gout(:,4)*180/pi,'m+','Markersize',3);
  	axis([X1 X2 -200 200]);		axis(axis);
     xlabel('������. �������, 1/�','Fontname','Arial Unicode MS'); 
     ylabel([num2str(kga) '*A, ��,  F, ����'],'Fontname','Arial Unicode MS'); grid;
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
      for j=1:2:(length(Vzon)-1)   % ������������� ������ 
         rab1=[Gout(Vzon(j),2) Gout(Vzon(j),2) Gout(Vzon(j+1),2) ...
                Gout(Vzon(j+1),2)];
         rab2=[130 150 150 130];
         semilogx(rab1,rab2,'r')
      semilogx(rab1,-rab2,'r')
      end
     	rab1=[int2str(Npp),'  ',int2str(Dat(1)),' '];
	rab2=[int2str(Dat(2)),' ',int2str(Dat(3))];
      %text(X1*(1+0.2*(log10(X2)-log10(X1))),-240,[rab1,rab2])
      text(X1*(X2/X1)^0.05,-240,[rab1,rab2]);
      title(tit)
	hold off
	axis('normal');
   else
   if Ngraf==2
   figure%(2)	
      Nf(2)=gcf;
      rab1=log10(nfc)-floor(log10(nfc));    % ���. �����. �� ��� �
      rab1=Vnorm(ceil(rab1*10)+1)*10^floor(log10(nfc)); X1=rab1;
      axis([0 rab1 -200 200]);	 axis(axis);
      Xmem(2)=0.8*X1;
	 % �������� ������
      plot(Gout(:,1),kga*Gout(:,3),'b',Gout(:,1),Gout(:,4)*180/pi,'m+',...
           Gout(:,1),100*log10(Gout(:,2)),'k-.','Markersize',3);
	axis([0 rab1 -200 200]);	 axis(axis);     
   xlabel('����� ����� N','Fontname','Arial Unicode MS'); 
   ylabel([num2str(kga) '*A, ��,  F, ����, 100*log10(om)'],'Fontname','Arial Unicode MS'); grid;
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
      for j=1:2:(length(Vzon)-1)   % ������������� ������ 
         rab1=[Gout(Vzon(j),1) Gout(Vzon(j),1) Gout(Vzon(j+1),1) ...
                Gout(Vzon(j+1),1)];
         rab2=[130 150 150 130];
         plot(rab1,rab2,'r')
         plot(rab1,-rab2,'r')
      end
      	rab1=[int2str(Npp),'  ',int2str(Dat(1)),' '];
	rab2=[int2str(Dat(2)),' ',int2str(Dat(3))];
      text(X1/10,-240,[rab1,rab2])
      title(tit)
	hold off
   axis('normal');
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure%(3)	
 Vaxis=[-90 40 -50 50];
      Nf(3)=gcf;
      for j=1:nfc
        if Gout(j,3)>-nidA		% Gout(j,3) - 20*log10(A)
           rab1=Gout(j,3)+nidA;
        else
           rab1=0;
        end
			Rgod(j)=rab1;
        XsmY(j,1)=rab1*cos(Gout(j,4))-nidA;
        XsmY(j,2)=rab1*sin(Gout(j,4));
      end
	

      rab1=pi/12;         % ��������� �����
      for j=1:25
         Ring(j,1:2)=[cos(j*rab1)  sin(j*rab1)];
      end
      rab1=nidA*ones(size(Ring(:,1)));
      rab2=(rab1+Rzon*Ring(:,1)).*cos(pi/6*Ring(:,2));
      rab3=(rab1+Rzon*Ring(:,1)).*sin(pi/6*Ring(:,2));
      axis(Vaxis);
	  axis(axis);
      hold on
  
      plot([Vaxis(1),Vaxis(2),Vaxis(2)],[Vaxis(4),Vaxis(4),Vaxis(3)],'k')   % �����
      patch(rab2-nidA,rab3,[0.9,0.9,0.9],'EdgeColor','none')	% ������ ����������� ����. �����
      patch((rab2-nidA)*2/3,rab3*2/3,'w','EdgeColor','none')	% ������ ����������� ����. �����

      rab1=[-80,-70,-50:10:-10,10,20,30]; rab2=zeros(size(rab1));		% ����� 10 �� +
      plot(rab1, rab2, '+k','Markersize',4)  %��������������
	  rab1=-20:10:20; rab2=-nidA*ones(size(rab1));	
      plot(rab2, rab1, '+k','Markersize',4)  %������������
      
      plot((Ring(1:23,1)-1)*nidA,Ring(1:23,2)*nidA,'xk','Markersize',5); % ������. ����. �������
      plot(0, 0, '+','Markersize',12);        plot(0, 0, 'o')   %����� � ����. �����

        %%%%%%%%%%%%%%%%%%%% �������� ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Mt1=-54: 12:-6; Mt2=43*ones(1,5);
      Mt=[Mt1,Mt1; Mt2,-Mt2];	% ������ ��������� ������� �������� ������
      Col=['r  ';'b: ';'k--';'m* ';'c--';'y: ';'g  '];
      Colar=['r';'b';'k';'m';'c';'y';'g'];	Ncl=7;
		s=0;   Vr=XsmY(1,:)';
      scol=rem(s,Ncl)+1;
for j=2:nfc
   rab1=~(((Rgod(j-1)<10)&(Rgod(j)>=10))|(j==nfc));
   %rab1=~(((Rgod(j-1)<10)&(Rgod(j)>=10))|...
		%((Rgod(j-1)<nidA)&(Rgod(j)>=nidA))|(j==nfc));
   if rab1	% ������������� -20�� �����
	   Vr=[Vr,XsmY(j,:)'];
	else
	   Vr=[Vr,XsmY(j,:)'];
	   scol=rem(s,Ncl)+1;
     %plot(Vr(1,:),Vr(2,:),'r'); 
     plot(Vr(1,:),Vr(2,:),Col(scol,:),'Markersize',3); % ��������� ����� �������� �����
	   s=s+1;
	   Vr=XsmY(j,:)';
   end
   if (Gout(j,2)-1)*(Gout(j-1,2)-1)<=0	%������� ����� om=1
      rab=((XsmY(j,1)-XsmY(j-1,1))^2+(XsmY(j,2)-XsmY(j-1,2))^2)^0.5;
       Ear=((XsmY(j,1)-XsmY(j-1,1))+i*(XsmY(j,2)-XsmY(j-1,2)))/rab;
            Rotarr=(arrow*Ear).';
       patch(XsmY(j,1)+real(Rotarr),XsmY(j,2)+imag(Rotarr),'r');	% ������� om=1
   end
	if j>2
      rab1=(Rgod(j-1)-Rgod(j-2))*(Rgod(j)-Rgod(j-1));
      if (rab1<0)&(Rgod(j)>10)	%�������� �������
         plot(XsmY(j,1),XsmY(j,2),'.')
            % �������:
            % Gout: [N,om,A,F(XsmY),Fsumangl(XsmY)]
            rab=abs(Gout(j,4)-Gout(j-1,4))<pi;
            rab1=-(~rab-rab)*pi/2;        
            if (Gout(j,4)<Gout(j-1,4))
               Ear=cos(Gout(j,4)-rab1)+i*sin(Gout(j,4)-rab1);
            else
               Ear=cos(Gout(j,4)+rab1)+i*sin(Gout(j,4)+rab1);
            end
            Rotarr=(arrow*Ear).';
            patch(XsmY(j,1)+real(Rotarr),XsmY(j,2)+imag(Rotarr),Colar(scol+1,:));	% �������	   
         if ~isempty(Mt)
            [L,k]=min((XsmY(j,1)-Mt(1,:)).^2+(XsmY(j,2)-Mt(2,:)).^2);
            plot([XsmY(j,1),Mt(1,k)],[XsmY(j,2),Mt(2,k)],':g')
            plot(Mt(1,k),Mt(2,k),'.')
            text(Mt(1,k),Mt(2,k)+3*sign(Mt(2,k)),num2str(round(Gout(j,2)*10)/10));
            rab=length(Mt(1,:));
         if rab==1
            Mt=[];
         else	%��������� �������������� �������
            switch k
            case rab,Mt=Mt(:,1:k-1);
            case 1,Mt=Mt(:,2:rab);
            otherwise,Mt=[Mt(:,1:k-1),Mt(:,k+1:rab)];
            end
         end
         end
      end
   end
 end

 
      text(25,-46,['Npp=',int2str(Npp)])
      if sum(size(XA))>0
         XAg=[round(XA(:,1)*10)/10 round(XA(:,2)*10)/10];
         text(14,43,'om        A')
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
      text(-80,-46,[int2str(Dat(1)) ' ' int2str(Dat(2)) ' ' int2str(Dat(3))])
      title(tit);
      hold off
      set(gca,'XTick',[],'YTick',[])
      axis('normal');
   end	%Ngraf==2	
end	%Ngraf==1
end	%for qgr=1:Numgraf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


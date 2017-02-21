function nichols_os(Name)
%Анализ систем с плотным спектром осцилляторов nichols_os(Name) 
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
disp(['Порядок системы ', Nnotes.Notes,' равен ',num2str(Nsys)])
if Ts>0
   Npp=sum(abs(Pol)>1)
else
   Npp=sum(real(Pol)>0)
end
%--------------------------------------------------------------------------%
figure
nichols(Name)
V=axis;
axis([V(1:2),-100,150])
title(tit);
hold on
text(0.85,0.96,['Npp=',int2str(Npp)],'units','normalized')
text(0.8,0.05,date,'units','normalized')
[mag,phase,w] = nichols(Name);
Nw=length(w);
Adb=20*log10(mag(:));
phase=phase(:);
%Разметка w
diap=log10([min(w),max(w)]);
rab=round(diap(1)*5)/5;
Sw(1)=rab;
s=1;    q=1;
while Sw(s)<=diap(2)
    if abs(round(Sw(s))-Sw(s))<0.1
        Sw(s)=round(Sw(s));
        Swc(q)=Sw(s);
        if abs(Swc(q))<0.1
            qe=q;
        end
        q=q+1;
    end
    s=s+1;
    Sw(s)=Sw(s-1)+0.2;
end
Sw=10.^Sw;
Swc=10.^Swc;
%Метки частот
SAdb=interp1(w,Adb,Sw);
Sphase=interp1(w,phase,Sw);
SAdbc=interp1(w,Adb,Swc);
Sphasec=interp1(w,phase,Swc);
plot(Sphase,SAdb,'.')
plot(Sphasec,SAdbc,'r.','Markersize',20)
text(Sphasec(qe)-(V(2)-V(1))/30,SAdbc(qe),'1')
%Определение запасов по амплитуде
s=1;    wA=[];
for j=2:Nw
    rab=round((phase(j)-180)/360)*360+180;
    if (phase(j)-rab)*(phase(j-1)-rab)<=0
        rab1=interp1([phase(j-1),phase(j)],[w(j-1),w(j)],rab);
        rab2=interp1([phase(j-1),phase(j)],[Adb(j-1),Adb(j)],rab);
        if rab2>=-30
            wA(:,s)=[rab1; rab2];
            s=s+1;
        end
    end
end
wA    
%Определение запасов по фазе
s=1;    wF=[];
for j=2:Nw
    rab=round((phase(j)-180)/360)*360+180;
    if Adb(j)*Adb(j-1)<=0
        rab1=interp1([Adb(j-1),Adb(j)],[w(j-1),w(j)],0);
        rab2=interp1([Adb(j-1),Adb(j)],[phase(j-1),phase(j)],0);
        if abs(rab2-rab)<90
            wF(:,s)=[rab1; rab2-rab];
            s=s+1;
        end
    end
end
wF    
% figure
% plot(w,phase)
% grid
%--------------------------------------------------------------------------%
if length(wA)>0
     text(0.05,0.96,'om','units','normalized')
     text(0.15,0.96,'A','units','normalized')
     for j=1:length(wA(1,:))
        text(0.05,0.95-j*0.05,num2str(round(wA(1,j)*10)/10),'units','normalized')
        text(0.15,0.95-j*0.05,num2str(round(wA(2,j)*10)/10),'units','normalized')
     end
end
if length(wF)>0
     text(0.3,0.96,'om','units','normalized')
     text(0.4,0.96,'F','units','normalized')
     for j=1:length(wF(1,:))
        text(0.3,0.95-j*0.05,num2str(round(wF(1,j)*10)/10),'units','normalized')
        text(0.4,0.95-j*0.05,int2str(wF(2,j)),'units','normalized')
     end
end

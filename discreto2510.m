clear 
close all
clc

% Parametros 
KLA = 4;       %%%%            
KNH=3; % (tambien=1) mg N1-1 coeficiente medio de saturacion de amonio para biomasa autrofa           
KNO = 0.5;     % mg N1-1 coeficiente medio de saturacion de nitrato para biomasa autotrofa      
KOA = 0.4;     % mg O21-1 coeficiente medio de saturacion de oxigeno para biomasa autotrofa  
KOH = 0.2;     % mg O21?1 coeficiente medio de saturacion de oxigeno para biomasa heterotrofa
KS = 20;       % mg COD1/1 Coeficiente medio de saturacion del sustrato para biomasa heterotrofa
KX = 0.03;     %(â€”)    
%condiciones iniciales
% SIin = 2;    % mg CODl/1 concentracion demateria organica soluble inerte       
SNDin = 9;     % mg Nl/1concentracion de material organico biodegradable soluble     
%SNDin = 15;
SNHin = 15;     % mg Nl?1 concentracion de nitrogeno de amonio soluble    
SNOin = 1.0;    % mg  Nl?1 concentracion de nitrato y nitrito nitrogenado    
SOin = 2.0;%%%%  % mg l?1 original 2 concentracion de oxigeno soluble       
SOmax = 10;  %%% % mg l?1 Concentracion maxima de oxigeno soluble     
SSin = 100;      % g CODl/l(200) concentracion de sustrato soluble facilmente biodegradable   
% SSin=200;          
XBAin = 0;        % mg CODl/l biomasa autotrofa activa 
XBHin = 0;  %%%   % mg COD l/l biomasa heterotrofa activa   
XNDin = 0;        % mg COD l/l Nitreno organico particulado biodegradable  
% XIin = 3.0;     % mg l/l Materia organica inerte particulada  
% XPin = 0;  

% % Biodegradables

XSin = 0.0;      %(en el otro tengo 100)   
% XSin=100;             
YA = 0.24;  %%%        
YH = 0.67;  %%%        
bA = 0.05;   %        
bH = 0.22;   %        
fp = 0.08;   %       
iXB = 0.086; %       
iXP = 0.06;  %        
kA = 0.081;  %       
kh = 3.0;    %       
ng = 0.8;    %        
nh = 0.4;    %       
mmaxA = 0.8; %%%       
mmaxH = 6;     %%%     
b = 0.24;    % day-1 (en el otro tengo 0.9)      

% q = 3000;          
% V = 1250;           
q=15; %2050;%7050;           
V = 15;    %7050;%2050;        %  Volumen del reactor
tau = V/q;       
% q=30; %caudal
% V=15;  %volumen    %

% Velocidades de reaccion correspondientes

% M2= Z(1)/(KS+Z(1));             % M2 = SS/(KS + SS);
% M8a = Z(5)/(KOA + Z(5));         % M8a =SO/(KOA + SO);
% M8h = Z(5)/(KOH + Z(5));        % M8a = SO/(KOH + SO);
% M9 = Z(6)/(KNO + Z(6));         % M9 = SNO/(KNO + SNO);
% M10 = Z(7)/(KNH + Z(7));         % M10 = SNH/(KNH + SNH);
% I8 = KOH/(KOH + Z(5));          % I8 = KOH/(KOH + SO);
% ksat = Z(2)/(KX*Z(3) + Z(2));   % ksat = XS/(KX*XBH + XS);
% qr = 0.5*q;                     % m3/dia

Ne=9; % numero de estados 
dt=.000001;% periodo de muestreo ( cada cuanto tomo muestras en el horizonte) tao
td = 0:dt:4; % vector de tiempo de simulacion (horizonte de tiempo)
ntotal=size(td,2); %para determinar la longitud del vector de tiempo 
% Z=zeros(Ne,ntotal); %para crear una matriz (para inicializar los estados )
% condiciones iniciales de los estados: 200 100 2250 250 0 0 15 9 0
Z=[200; 100; 2250; 250; 0; 0; 15; 9; 0];
% Z= [NE1,ntotal]; 
%Z(2,1)XS
%Z(3,1)XHB
%Z(4,1)XBA
%Z(5,1)SO
%Z(6,1)SNO
%Z(7,1)SNH
%Z(8,1)SND
%Z(9,1)XND


for k=1:ntotal-1
    
M2= Z(1,k)/(KS+Z(1,k));            
M8a = Z(5,k)/(KOA + Z(5,k));         
M8h = Z(5,k)/(KOH + Z(5,k));        
M9 = Z(6,k)/(KNO + Z(6,k));         
M10 = Z(7,k)/(KNH + Z(7,k));         
I8 = KOH/(KOH + Z(5,k));          
ksat = Z(2,k)/(KX*Z(3,k) + Z(2,k));   
qr = 0.5*q;                    

    

    
    dz1=((q/V)*(SSin - Z(1,k)))-((mmaxH/YH)*M2*(M8h+I8*M9*ng)*Z(3,k))+(kh*ksat*(M8h+nh*I8*M9)*Z(3,k));
    dz2=((q/V)*(XSin - Z(2,k)))+ ((qr/V)*(b-1)*Z(2,k))+((1-fp)*(bH*Z(3,k)+bA*Z(4,k)))-(kh*ksat*(M8h+nh*I8*M9))*Z(3,k);
    dz3=((q/V)*(XBHin-Z(3,k)))+((qr/V)*(b-1)*Z(3,k))+ (mmaxH*M2*M8h*Z(3,k))+ (mmaxH*M2*I8*M9*ng*Z(3,k))-(bH*Z(3,k));
    dz4=((q/V)*(XBAin-Z(4,k)))+ ((qr/V)*(b-1)*Z(4,k))+ (mmaxA*M10*M8a*Z(4,k))-(bA*Z(4,k));
    dz5=((q/V)*(SOin-Z(5,k)))+ (KLA*(SOmax-Z(5,k)))-(((1-YH)/YH)*mmaxH*M2*M8h*Z(3,k))-(((4.57-YA)/YA)*mmaxA*M10*M8a*Z(4,k));
    dz6=(q/V)*(SNOin-Z(6,k))-((1-YH)/(2.86*YH))*mmaxH*M2*I8*M9*ng*Z(3,k)+(1/YA)*mmaxA*M10*M8a*Z(4,k);
    dz7=(q/V)*(SNHin-Z(7,k))-iXB*mmaxH*M2*(M8h+I8*M9*ng)*Z(3,k)-(iXB+(1/YA))*mmaxA*M10*M8a*Z(4,k)+kA*Z(8,k)*Z(3,k);
    dz8=(q/V)*(SNDin-Z(8,k))- kA*Z(8,k)*Z(3,k)+kh*ksat*(M8h+nh*I8*M9)*Z(3,k)*(Z(8,k)/Z(2,k));
    dz9=(q/V)*(XNDin-Z(9,k))+(qr/V)*(b-1)*Z(9,k)+(iXB-fp*iXP)*(bH*Z(3,k)+bA*Z(4,k)) - kh*ksat*(M8h+nh*I8*M9)*Z(3,k)*(Z(9,k)/Z(2,k));
     Z(:,k+1)=Z(:,k)+dt*([dz1;dz2;dz3;dz4;dz5;dz6;dz7;dz8;dz9]);%utilizando euler hacia adelante 
     
     
end
figure (1)
% subplot(3,3,1)
plot(td,Z(1,:),'.k'),xlabel('tiempo'),ylabel('(SS)Sustrato soluble(mg/L)')
% grid on
hold on 

figure (2)
% subplot(3,3,2)
plot(td,Z(2,:)),xlabel('tiempo'),ylabel('(XS)Sustrato particulado lentamente biodegradable')
% grid on
hold on

figure (3)
% subplot(3,3,3)
plot(td,Z(3,:)),xlabel('tiempo'),ylabel('(XBH) Biomasa activa particulada Heterotrofa')
% grid on
hold on

figure (4)
% subplot(3,3,4)
plot(td,Z(4,:)),xlabel('tiempo'),ylabel('(XBA) Biomasa activa particulada heterotrofa')
% grid on
hold on

figure (5)
% subplot(3,3,5)
plot(td,Z(5,:)),xlabel('tiempo'),ylabel('(SO) Oxigeno soluble')
% grid on
hold on

figure (6)
% subplot(3,3,6)
plot(td,Z(6,:)),xlabel('tiempo'),ylabel('(SNO)Nitrato soluble y nitrito')
% grid on
hold on

figure (7)
% subplot(3,3,7)
plot(td,Z(7,:)),xlabel('tiempo'),ylabel('(NH4 Y NH)Amonio soluble y nitrogeno')
% grid on
hold on

figure (8)
% subplot(3,3,8)
plot(td,Z(8,:)),xlabel('tiempo'),ylabel('(SND)Nitrogeno organico soluble biodegradable')
% grid on
hold on

figure (9)
% subplot(3,3,9)
plot(td,Z(9,:)),xlabel('tiempo'),ylabel('(XND)Nitrogeno Organico particulado biodegradable')
% grid on
hold on

figure (10)
% subplot(3,3,)
plot (td,Z(1,:)+Z(2,:)+5);xlabel('tiempo'),ylabel('(COD Total)Dinamica de la demanda quimica de oxigeno')
% grid on 






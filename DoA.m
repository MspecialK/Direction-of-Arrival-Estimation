%% RADAR IMAGING
%  Mykhaylo Kotsaba

close all;
clear all;
clc;

%  TASK 2 - Implementation of DOA estimation by 2D-FFT
%% A.                          
f0=9.6e9;
lambda = 3e8/f0; 
ds=lambda/5;

%% B.
limit=[10, 20]; 
p=[0,0,0];                  

% We create a source point in our X,Y,Z limits coordinations.
while ((-limit(1)<=p(1) && p(1)<limit(1) ) || (-limit(1)<=p(2) && p(2)<limit(1)) || (-limit(1)<=p(3) && p(3)<limit(1) )) 
    p=[2*limit(2)*rand-limit(2), 2*limit(2)*rand-limit(2), limit(2)*rand];           
end
real_distance=norm(p);                                                    
real_phi=acos(p(3)/real_distance);                                        
real_teta=atan(p(1)/p(2)); 
Lambdax=lambda/sin(real_phi)/cos(real_teta);
Lambday=lambda/sin(real_phi)/sin(real_teta);

N=max(round(Lambdax/ds),round(Lambday/ds));
[X,Y]=meshgrid(((-(N-1)/2):(N-1)/2)*ds,((-(N-1)/2):(N-1)/2)*ds);


% 3D Plots of the configuration
figure;
subplot(1,2,1)   
hold on;
title('Position of antennas and of the source');
plot3(X,Y,zeros(size(X)),'.r');
plot3(p(1),p(2),p(3),'*b'); 
axis([-limit(2) limit(2) -limit(2) limit(2) 0 limit(2)]);
xlabel('x-axis [m]');ylabel('y-axis [m]');zlabel('z-axis [m]');

% Drawing some lines to make the graphic more clear
line([p(1),0],[p(2),0],[p(3),0],'LineWidth',0.2,'Color','k');        % distanza
line([p(1),p(1)],[p(2),p(2)],[p(3),0],'LineWidth',0.2,'Color','b');  % distanza asse z 
line([p(1),0],[p(2),p(2)],[0,0],'LineWidth',0.2,'Color','r');        % distanza asse x
line([p(1),p(1)],[p(2),0],[0,0],'LineWidth',0.2,'Color','g');        % distanza asse y
line([-limit(2),limit(2)],[0,0],[0,0],'LineWidth',0.2,'Color','k');  % asse x
line([0,0],[-limit(2),limit(2)],[0,0],'LineWidth',0.2,'Color','k');  % asse y
hold off;
grid on;

% Second graphic to show in detail the array
subplot(1,2,2)
plot3(X,Y,zeros(size(X)),'*r');
title('Array of NxN Anthennas');
xlabel('x-axis [m]');ylabel('y-axis [m]');zlabel('z-axis [m]');                             
grid on;


%% C. Received Signal

Vamp=1;                                                                    % Signal emitted amplitude
distance_m_n=sqrt((p(1)-X).^2 +(p(2)-Y).^2 +(p(3)).^2);                    % Distance between the point P and the (n,m) antenna
TSignal = Vamp./distance_m_n.*exp(-1j.*(2*pi/lambda).*distance_m_n);       % Signal received. 

% Plot of the received Signal Time Domain
figure
surf(X,Y,real(TSignal)); 
xlabel('x-axis [m]');ylabel('y-axis [m]');zlabel('Signal received ');
title('Signal received on the array matrix');

%% D. DFT Analysis

Ns = 128;
df=(1/ds)/Ns; 
FSignal=fftshift(fft2(TSignal,Ns,Ns)); 
[frx, fry]=meshgrid(((-(Ns)/2):(Ns-1)/2)*df, ((-(Ns)/2):(Ns-1)/2)*df); 

MaxVal = max(max(abs(FSignal)));                      % Max value of the dft
[xm,ym]=find(abs(FSignal)==MaxVal);                   % Index of the max value
Fx=frx(1,xm);                                         % Fx value 
Fy=fry(ym,1);                                         % Fy value

% Plot Signal in Frequency Domain
figure
surf(frx,fry,abs(FSignal)/N/N);xlabel('x [1/m]');ylabel('y [1/m]');zlabel('Sx');title('Fourier 2D Traansform on the received Signal')
teta=atan(Fy/Fx);
phi=asin(Fx*lambda/cos(teta));

%% E. NOISE Analysis

SNR=100;  

Amp=sqrt(1/SNR)*rand(N);  
phase=2*pi*rand(N);
TSignalNoise=TSignal+Amp.*exp(-1j.*phase);
FSignalNoise=fftshift(fft2(TSignalNoise,Ns,Ns))/N/N;      

MaxValNoise = max(max(abs(FSignalNoise)));                     % Max value of the dft with noise
[xmnoise,ymnoise]=find(abs(FSignalNoise)==MaxValNoise);        % Index of the max value noise
FxNoise=frx(1,xmnoise);                                        % Fx value noise
FyNoise=fry(ymnoise,1);                                        % Fy value Noise

tetaNoise=atan(FyNoise/FxNoise);
phiNoise=asin(FxNoise*lambda/cos(tetaNoise));

% Plotting Time and Frequency Domain
figure
subplot(1,2,1);
surf(X,Y,real(TSignalNoise)); 
xlabel('x-axis [m]');ylabel('y-axis [m]');zlabel('Signal received ');
title('Signal received on the array matrix with noise');
subplot(1,2,2);
surf(frx,fry,abs(FSignalNoise));
xlabel('x [1/m]');ylabel('y [1/m]');zlabel('Sx');title('Fourier 2D Traansform on the received Signal + Noise')



%% F.
%  Second target 
p2=[0,0,0];

while ((-limit(1)<=p2(1) && p2(1)<limit(1) ) || (-limit(1)<=p2(2) && p2(2)<limit(1)) || (-limit(1)<=p2(3) && p2(3)<limit(1) )) 
    p2=[2*limit(2)*rand-limit(2), 2*limit(2)*rand-limit(2), limit(2)*rand];                   
end

real_distance2=norm(p2);                                                        % Real distance from the [0,0,0]
real_phi2=acos(p2(3)/real_distance2);                                           % Real elevation angle of the point
real_teta2=atan(p2(2)/p2(1)); 

% Signal of 2 sources
distance_m_n2 = sqrt((p2(1)-X).^2 +(p2(2)-Y).^2 +(p2(3)).^2);    
TSignal2 = Vamp./distance_m_n2.*exp(-1j.*(2*pi/lambda).*distance_m_n2)+Vamp./distance_m_n.*exp(-1j.*(2*pi/lambda).*distance_m_n);    % Signal received. the amplitude is negligible for our purpose setted at 1 divided per distance

Distancepp=sqrt( (p(1)-p2(1)).^2+(p(2)-p2(2)).^2+(p(3)-p2(3)).^2);

%Plot
figure
FSignal2=fftshift(fft2(TSignal2,Ns,Ns))/N/N; 
surf(frx,fry,abs(FSignal2));


%% TASK 3 
% A. 

dx=0.5;
space_dx=-limit(2):dx:limit(2);
dim=size(space_dx,2);

[XX,YY,ZZ]=meshgrid(space_dx);
BackSignal=zeros(dim, dim, dim);

%% B.

for ii=1:dim
    for jj=1:dim
        for hh=1:dim
            DistPositive=sqrt((X-space_dx(ii)).^2 +(Y-space_dx(jj)).^2 +(space_dx(hh)).^2);    
            BackSignal(jj,ii,hh)=sum(TSignal.*DistPositive.*exp(1j*(2*pi/lambda)*DistPositive),'all'); 
        end
    end
end

xslice=[p(1)];  %% Slicing in the real p point 
yslice=[p(2)];
zslice=[p(3)];

figure
slice(XX,YY,ZZ,abs(BackSignal),xslice,yslice,zslice);
xlabel('x-axis [m]');ylabel('y-axis [m]');zlabel('Signal received ');




fprintf("Real paramaeters:\n");
fprintf("Phi: %.2f gradi \nTeta %.2f gradi\n----------------------\n",abs(rad2deg(real_phi)),rad2deg(real_teta) );
fprintf("Calculated :\nPhi: %.2f gradi \nTeta %.2f gradi\n----------------------\n", abs(rad2deg(phi)), rad2deg(teta));
fprintf("Calculated with Noise:\nPhi: %.2f gradi \nTeta %.2f gradi\n----------------------\n", abs(rad2deg(phiNoise)), rad2deg(tetaNoise));




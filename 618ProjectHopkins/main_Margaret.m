%main file for MUMT 618 Final project, Margaret Hopkins
%this function is a modification of an implementation of Adachi and Sato (1996), written by Vincent Freour
%last modified December 13, 2022 by Margaret Hopkins

function [t, p]=main_Margaret(Eeqy, flip, pmax);
%%%%%%%%SETUP%%%%%%%%%%%%
%clearvars -except g Eeqy

dir=genpath('.');
addpath(dir)

%pmax = 3500; %in pascals, change for the different partials
N = 24000; %number of samples to compute (2s duration)
fs = 8000; %sampling frequency


a = 'Impedance_gap1.txt';
a = importdata(a);
f = a(:,1);
ReZ = a(:,2);
ImZ = a(:,3);
Z = ReZ + 1i*ImZ;
r=trumpet_imp(Z, f,fs,N);


clearvars -except pmax r phase_offset N fs Eeqy pmax flip

%%input parameters

flip=ones(1,N)*flip; %*2^(1./12.);     % lip resonance frequency, vector of sample values (60-700Hz)
%can multiply by  *2^(1./12.) to get other notes

%lip parameters
Q = 3.0; %from 1996 paper, quality factor
b = 7E-3; %lip breadth
d = 2E-3; %thickness of lips
%geometrical parameters

%lip joint position from 1996: "E" corresponds to xi in the paper
Ejointx = 0; %x coordinate of B (constant position of "hinge" at top/bottom of lip)
Ejointy = 4E-3; %y coordinate of B
Eeqx = 1E-3; %x coordinate of lip rest position (C)
%Eeqy = 0.3E-3; %y coordinate of lip rest position (C) -0.1E-3 to -2.0E-3

%parameters
T = 1/fs; %sample period
rho = 1.1769; %density of air (kg/m^3)
c = 347.23; %speed of sound in air (m/s) at 27.03 degrees C

%mouthpiece geometry
diam = 0.0248;
%diam = 0.01685; %diameter of trumpet mouthpiece - I forget which one (Bach 5C?)
Scup = pi*(diam/2)^2; %area of mouthpiece entryway

%lip mechanical properties
%k is proportional to flip because flip = sqrt(k/m)/(2pi)
m = 1.5./(flip.*(2*pi).^2); %transverse model in 1995, whole model in 1996
k = 1.5*flip;

%characteristic impedance
Zc = rho*c/Scup;

%time vector
%t = 0:T:N*T-T;

%initial conditions - start at equilibrium for first two terms
Ex(1) = Eeqx;
Ex(2) = Eeqx;
Ex = [Ex zeros(1, N-2)];

Ey(1) = Eeqy;
Ey(2) = Eeqy;
Ey = [Ey zeros(1, N-2)];

Slip(1) = max([2*b*Ey(1) 0]); %initial cross-sectional area of lip orfice, 0 if the opening is "negative"
Slip(2) = max([2*b*Ey(2) 0]); 
Slip = [Slip zeros(1, N-2)]; %initializing Slip with initial values and zeros

%initializing variables with zeros
v=zeros(1,N);       % velocity in x direction
z=zeros(1,N);       % velocity in y direction
Uac=zeros(1,N);     % time varying component of volume velocity
Ulip=zeros(1,N);    % volume velocity at lip opening
U=zeros(1,N);       % volume velocity?
p=zeros(1,N);       % pressure
plip=zeros(1,N);    % average pressure at lip opening

%blowing pressure
%p0 = [pressure_input(pmax,N,0.1,0.1,fs)];
p0 = [pressure_input(pmax,N,(500*fs/48000)/fs,(500*fs/48000)/fs,fs)];


%ramp up to pmax for 500 samples, ramp down for 500 samples at the end

%[amp, ~] = onset_tuning_1(0, 1, 0.1, 1, (N-Nz)/fs, fs);
%onset tuning stuff that doesn't get used(?)

C = 0;
index = 1;
k2 = 1;
for n = 1:N-1

    %lip mechanical properties
    wL(n) = 2*pi*flip(n); %lip frequency into angular speed
    mL(n) = 1.5/((2*pi)^2*flip(n)); %lip mass in kg from table 1 in 1996 paper
    kL(n) = 1.5*flip(n); %stiffness of lips, in N/m

    M(n)=sqrt(mL(n)*kL(n))/Q; % coefficient on dXi/dt without 1/2 factor in 
    % eq 1 1996

    %%%%%%%%Step 1 in 1996 paper: finding new xi (E) at one step ahead%%%%%

    Ex(n+1) = T*v(n)+Ex(n); %adding how much Ex (lip x position Xi_x) has 
    % moved since last timestep based on previous x velocity

    v(n+1)=(2*T/mL(n))*((mL(n)/2/T - 0.5*M(n))*v(n) - 0.5*kL(n)*(Ex(n)-Eeqx) + b*(p0(n)-p(n))*(Ejointy-Ey(n)));
    %next x-velocity value: terms are very similar to eq (1) in 1996, but
    %missing F_bernoulli, with extra term of (2T/m)*(m/2T), and with just
    %x-direction in F_restore term, and just y-direction length in F_delta
    %term, and x-velocity instead of dXi/dt

    Ey(n+1)=T*z(n)+Ey(n); %amount y moves to next term (z is y-velocity)
    
    z(n+1)=(2*T/mL(n))*((mL(n)/2/T - 0.5*M(n))*z(n) - 0.5*kL(n)*(Ey(n)-Eeqy) + b*(p0(n)-p(n))*Ex(n) + b*d*plip(n) - C*kL(n)*(Ey(n)-Eeqy));
    %next y-velocity value: terms the same as above for v(n+1), but with
    %y's and x's swapped and an extra term that is set to 0 since C=0

    %%%%%%%%%%%%Step 2 in paper: Calculating a new S_lip and U_lip%%%%%%%%%

    Slip(n+1)= max([2*b*Ey(n+1) 0]); %from eq5 in 1996: when lips are closed it is set to 0



    %if Slip(n+1)<1E-16
        %Q=0.5;
    %else
        %Q=3;
    %end

    if Slip(n+1) == 0 %i.e. lip area == 0 
        Ey(n+1) = 0; % set the y value to 0 (y "height" of lips is 0) so that position does not become negative
        Ulip(n+1) = b*((Ex(n)-Ejointx)*((Ey(n+1)-Ey(n))/T)-(Ey(n)-Ejointy)*((Ex(n+1)-Ex(n))/T)); %volume velocity from (6) in 1996 paper with discrete approximation of derivatives
    else
        Ulip(n+1) = b*((Ex(n)-Ejointx)*((Ey(n+1)-Ey(n))/T)-(Ey(n)-Ejointy)*((Ex(n+1)-Ex(n))/T)); %volume velocity from (6) in 1996 paper with discrete approximation of derivatives


    %%%Step 3 in paper:Calculate U_acoust and p from (7), (8), and (9)%%%%%
        if sign(p0(n)-p(n))== 1
             Uac(n+1) = Slip(n+1)*sign(p0(n)-p(n))*sqrt((2*(p0(n)-p(n)))/rho);
        
        elseif sign(p0(n)-p(n))==0
           Uac(n+1) = Slip(n+1)*sign(p0(n)-p(n))*sqrt((2*(p0(n)-p(n)))/rho);
       
        elseif sign(p0(n)-p(n))==-1
           Uac(n+1) = Slip(n+1)*sign(p0(n)-p(n))*sqrt((2*(p(n)-p0(n)))/rho);
           
        end %acoustic volume velocity
    %seems to be from equations 7 and 8 in the paper, not sure exactly how
    %Freour's code has conditions on the sign of p0(n)-p(n), but the result
    %is the same for all 3 set of conditions
    end

    %summing lip and acoustic volume flow to get total

    U(n+1) = Ulip(n+1)+Uac(n+1);

    %equation 9 in the paper: this is where measured impedance comes in

    h=1;
    for i=0:n-1
        Y(h) = (Zc*U(n-i)+p(n-i))*r(i+1); 
        %r is the reflectance from the measured impedance
        %Y(n) is the quantity that is integrated in equation (9)
        h=h+1;
    end

    p(n+1) = trapz(Y) + Zc*U(n+1); %calculation of (9)
    %trapz is numerical integration
end

%%plotting
% gray = [0.6 0.6 0.6];
% 
% figure(1)
% plot(p0-p)
% title('Delta p')
% hold off
% 
% figure(2)
% plot(p)
% hold on
% legend('p');
% hold off
% 
% figure(3)
% plot(Ex, Ey)
% xlim([0 3.5E-3])
% ylim([0 3.5E-3])
% title('Ex vs. Ey')
% 
 t=0:T:N*T-T;
% 
% t1 = 0.18;
% t2=0.23;
% 
% t3=0.7;
% t4=0.75;
% 
% t5=1.4;
% t6=1.45;
% 
% t7=1.8;
% t8=1.85;
% hold off
% 
% figure(4)
% subplot(2,4,1:4)
% plot(t,p,'k')
% hold on
% xlabel('Time (s)','fontsize',14)
% ylabel('Pa','fontsize',14)
% hold on
% line([t1 t1],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t2 t2],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t3 t3],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t4 t4],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t5 t5],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t6 t6],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t7 t7],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% hold on
% line([t8 t8],[-1.5E4 1E4],'color','k','linestyle','--','linewidth',1,'color','r')
% ylim([-1.5E4 1E4])
% subplot(2,4,5)
% plot(Ex(t1*fs:t2*fs),Ey(t1*fs:t2*fs),'k')
% xlabel('\xi_x (m)','fontsize',14)
% ylabel('\xi_y (m)','fontsize',14)
% ylim([0 2.5E-3])
% xlim([-0.003 0.008])
% subplot(2,4,6)
% plot(Ex(t3*fs:t4*fs),Ey(t3*fs:t4*fs),'k')
% xlabel('\xi_x (m)','fontsize',14)
% ylabel('\xi_y (m)','fontsize',14)
% ylim([0 2.5E-3])
% xlim([-0.003 0.008])
% subplot(2,4,7)
% plot(Ex(t5*fs:t6*fs),Ey(t5*fs:t6*fs),'k')
% xlabel('\xi_x (m)','fontsize',14)
% ylabel('\xi_y (m)','fontsize',14)
% ylim([0 2.5E-3])
% xlim([-0.003 0.008])
% subplot(2,4,8)
% plot(Ex(t7*fs:t8*fs),Ey(t7*fs:t8*fs),'k')
% xlabel('\xi_x (m)','fontsize',14)
% ylabel('\xi_y (m)','fontsize',14)
% ylim([0 2.5E-3])
% xlim([-0.003 0.008])
% hold off
% 
% figure(6)
% t=(1:size(r))*1000./fs;
% plot(t,r)
% axis([0 100 -0.02 0.01])
% xlabel('Time (ms)','fontsize',14)
% ylabel('Refl. func.','fontsize',14)
% 

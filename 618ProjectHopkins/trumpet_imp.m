function [r]=trumpet_imp(Z, f, fs, N);

%load('C:\Users\Vincent\Documents\ISMA_2010\VT_impedance\Trombone\Trombone t-slide 0.mat');
%load('C:\Users\Vincent\Documents\ISMA_2010\VT_impedance\Trombone\Trombone_slide_0.mat');
%outputs: r is time-domain reflection function
%ir is ifft of I2
%I2 is impedance after smoothing
%ft sampled time indices
%R2 is reflection coefficient (after symmetrization?)
% 
% fs = 12000;
% N = 24000;
% a = 'Impedance_gap1.txt';
% a = importdata(a);
% f = a(:,1);
% ReZ = a(:,2);
% ImZ = a(:,3);
% Z = ReZ + 1i*ImZ;


        Z=double(Z); %convert to double class
        %length correction due to diameter change between probe and mouthpiece
        %Z = addpipe( Z, f, -0.0015 );

        %Z=Z*((0.01685/2)/0.00635)^2; %Diameter variation (probe to mouthpiece) correction ** modified for trumpet: uses radius

        %length correction due to the position of the lips downward the
        %mouthpiece - ** is this delta_x? where to find this?
        I = addpipe( Z, f, -0.003 );
        
        
%         rho = 1.1769;           % density of air (kg/m^3)
%         c = 347.23;             % speed of sound in air (m/s)
%         %Mouthpiece geometry
%         diam=0.0248;
%         Scup= pi*(diam/2)^2;
%         Zc=rho*c/Scup;
%         I=I*Zc;
        
        %I(end+1)=0;

       %f = [((0:249)*(f(2)-f(1)))'; f]; %for plotting fitted I below
       %against frequency

        %is this to get rid of noise in the early part of the impedance?
        %linear fitting of real part of I (impedance?) from 1 up to 58 -
        %for trumpet this just smooths a tiny bit at the bottom
        thres=400;
        reg=polyfit([1 thres],[0 real(I(thres-250))],1);
        realI(1:thres)=[1:thres]*reg(1)+reg(2);
        
        %linear fitting for imaginary part of I
        reg=polyfit([1 thres],[0 imag(I(thres-250))],1);
        imagI(1:thres)=[1:thres]*reg(1)+reg(2);
        
        I = [realI(1:250)'+j*imagI(1:250)'; I];
        I(251:thres) = realI(251:thres)+1i*imagI(251:thres);
       % I(1:thres)=realI(1:thres)+j*imagI(1:thres);
       %[pks, locs] = findpeaks(abs(I), 'MinPeakProminence', 20) %finding impedance
       %peaks
       %f(locs) gives frequencies of peaks
        
        R=((I-1)./(I+1)); %reflection function
   
   %window for fft
   %Half sine windowing of the high frequency part of R
   %d=5000;
%    d=1500*2;
%    u=0.001;
%    W=ones(1,length(R));
%    k=u;
%         for g=d:floor(d+pi/u);
%             W(g)=W(d)*((sin(pi/2-k)+1)/2);
%             k=k+u;
%         end
%         for g=floor(d+pi/u+1):length(R);
%             W(g)=0;
%         end
%    
%    R=R.*W';     
        
  %Symetrisation of R before inverse fft
fN=fs/2;
step= 24000/length(R);
R=R(1:floor(fN/step));

R2 = [R ;flipud(conj(R(2:end-1,:)))];
R2=resample(R2,N,length(R2));
r = real(ifft(R2));
%r(1)=0;
%%

%Causality correction
%tells program which side of real part the imaginary part is on
rc=ifftshift(r);
rc_neg=rc(1:floor(length(rc)/2));% added the floor on 10/04/2013
rc_pos=rc(floor(length(rc))/2+1:end);

for n=1:(length(rc)/2); 
rc_pos(n)=rc_pos(n)+rc_neg(floor(length(rc)/2)-n+1); 
rc_neg(n)=0;
end

r=rc_pos;

%ft=0:24000/length(I):24000-24000/length(I);

%Resampling

%r=resample(r,fs,48000); %uncomment if no resampling applied to R2
% r_s=resample(r,fs,48000);
% r=[r_s ;zeros(N-length(r_s),1)];

%setting endpoints to 0 (dc component)
r(1)=0;
r(end)=0;

%%
% commented that part on 10/04/2013

   %Half sine windowing of the high frequency part of r
%    d=0.2*fs;
%    u=0.001*48000/fs;%Check this...
%    W=ones(1,length(r));
%    k=u;
%         for g=d:floor(d+pi/u);
%             W(g)=W(d)*((sin(pi/2-k)+1)/2);
%             k=k+u;
%         end
%         for g=floor(d+pi/u+1):length(r);
%             W(g)=0;
%         end
%    display(size(W))     
%    r=r.*W';  

% r(0.1*fs:length(r))=0;

if N > length(r);
    r=[r ; zeros(N-length(r),1)];
    elseif N<length(r);
    r=r(1:N);
end


I2 = [I ;flipud(conj(I(2:end-1,:)))];
%I2=resample(I2,fs,48000);
ir = real(ifft(I2));

ft=0:48000/length(I2):48000-48000/length(I2);
        
% ir=resample(ir,fs,48000);
% 
% I2=resample(I2,48000,length(I2));
% 
% ft=0:length(I2)/48000:48000-length(I2)/48000;

%plot(abs(I2))


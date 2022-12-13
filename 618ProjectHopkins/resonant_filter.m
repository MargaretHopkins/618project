function [y]=resonant_filter(x,r,fc,fs,option);

a(1)=1;
a(2)=-2*r*cos(2*pi*fc/fs);
a(3)=r^2;

b(1)=(1-r^2)/2;
b(2)=0;
b(3)=-b(1);

if option=='filter';
    y=filter(b,a,x);
else
    
y=filtfilt(b,a,x);
end
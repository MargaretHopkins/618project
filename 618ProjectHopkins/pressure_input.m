function p=pressure_input(pmax,N,onset,decay,fs)

%d = onset duration
%N = total length

%fs = 48000;
%pmax = 6000;
%N = 96000;
%onset = (500*fs/48000)/fs;
%decay = (500*fs/48000)/fs;

p=zeros(1,N);

od=floor(onset*fs);
dd=floor(decay*fs);
sd=N-od-dd;

for n=1:od+1;
    
    p(n)=pmax*((sin((n-1)*pi/od-pi/2)+1)/2);

end

for n=od+2:od+sd;
    
    p(n)=pmax;
    
end


for n=od+sd+1:N;
    
    p(n)=pmax*((sin((N-n)*pi/dd-pi/2)+1)/2);

end

%plot(p)
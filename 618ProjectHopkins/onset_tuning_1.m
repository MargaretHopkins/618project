function [A,t]=onset_tuning_1(init,max,Q,steady,duration,fs);

N1=fs*Q/2;

t=0:1/fs:duration-1/fs;

for n=1:N1+1;
    
    A(n)=(max-init)*((sin((n-1)*pi/N1-pi/2)+1)/2)+init;
    
end

for n=N1+2:Q*fs;
    
    A(n)=(max-steady)*((sin((n-1)*pi/N1-pi/2)+1)/2)+steady; %May have to be corrected (cf matlab file pressure_input.m

end

for n=Q*fs+1:duration*fs;
    A(n)=steady;
end


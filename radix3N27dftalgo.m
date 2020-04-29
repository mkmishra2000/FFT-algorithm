%Radix 3 fft algo for N =27
function SEQ = radix3N27dftalgo(seq)
    SEQ = zeros(1,length(seq));
    
    %divide the sequence into three sequences
    %of length 9.
    x9f=zeros(1,length(seq)/3); %first 9 pt sequence
    x9s=zeros(1,length(seq)/3); %second 9 pt sequence
    x9t=zeros(1,length(seq)/3); %third 9 pt sequence

    %divide the sequence and store into another array
    %using decimation in time domain
    for i=1:9
        x9f(i)= seq(3*i-2);
        x9s(i)= seq(3*i-1);
        x9t(i)= seq(3*i);
    end
    
    %vacant arrays for storing the 9 pt dft
    X9f = zeros(1,length(seq)/3);
    X9s = zeros(1,length(seq)/3);
    X9t = zeros(1,length(seq)/3);
    
    X9f = radix3N9dftalgo(x9f);
    X9s = radix3N9dftalgo(x9s);
    X9t = radix3N9dftalgo(x9t);
    
    W27 = exp((-1j*2*pi)/27);
    
    for p =1: length(seq)/3
        SEQ(p)   = X9f(p) + W27^(p-1) * X9s(p) + W27^(2*(p-1)) * X9t(p);
        SEQ(p+9) = X9f(p) + W27^(p+8) * X9s(p) + W27^(2*(p+8)) * X9t(p);
        SEQ(p+18)= X9f(p) + W27^(p+17) * X9s(p) + W27^(2*(p+17)) * X9t(p);
    end
end
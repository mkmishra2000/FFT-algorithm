%9 pt DFT using 3pt DFT
function SEQ=radix3N9dftalgo(seq)
    SEQ = zeros(1,length(seq));
    
    x3f=zeros(1,length(seq)/3); %first 3 pt sequence.
    x3s=zeros(1,length(seq)/3); %second 3 pt sequence.
    x3t=zeros(1,length(seq)/3); %third 3 pt sequence.

    %distribution of sequence's element. 
    for i=1:length(seq)/3
        x3f(i) = seq(3*i-2);
        x3s(i) = seq(3*i-1);
        x3t(i) = seq(3*i);
    end

    %Vacant arrays for 3 pt DFT 
    Xf=zeros(1,length(seq)/3);
    Xs=zeros(1,length(seq)/3);
    Xt=zeros(1,length(seq)/3);
    
    %3pt dft
    Xf = esd113ptdtf(x3f);
    Xs = esd113ptdtf(x3s);
    Xt = esd113ptdtf(x3t);
    
    %9 pt twiddle factor
    W9 = exp((-1j*2*pi*1)/9);

    for p =1: length(seq)/3
        SEQ(p) = Xf(p) + W9^(p-1) * Xs(p) + W9^(2*(p-1)) * Xt(p);
        SEQ(p+3)=Xf(p) + W9^(p+2) * Xs(p) + W9^(2*(p+2)) * Xt(p);
        SEQ(p+6)=Xf(p) + W9^(p+5) * Xs(p) + W9^(2*(p+5)) * Xt(p);
    end

end
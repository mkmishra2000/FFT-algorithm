%3pt DFT of the sequence
function SEQ = esd113ptdtf(seq)
    SEQ = zeros(1,length(seq));
    
    W = exp((-1j*2*pi)/3);
    
    for i =1:length(seq)
        SEQ(i) = seq(1) + W^(i-1)*seq(2) + W^(2*(i-1))*seq(3);
    end 
end
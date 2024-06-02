function W = mywigner(Ex)
    N = length(Ex);														
    x = ifftshift(((0:N-1)'-N/2)*2*pi/(N-1));							
    X = (0:N-1)-N/2;
    EX1 = ifft( (fft(Ex)*ones(1,N)).*exp( 1i*x*X/2 ));					
    EX2 = ifft( (fft(Ex)*ones(1,N)).*exp( -1i*x*X/2 ));					
    W = real(fftshift(fft(fftshift(EX1.*conj(EX2), 2), [], 2), 2));		
end
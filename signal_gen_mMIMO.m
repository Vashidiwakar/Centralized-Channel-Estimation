% Signal Generation
function Yk = signal_gen_mMIMO(Phi, Xk, M, T, sigma2_k, A_R) 
    Nk = sqrt(sigma2_k/2)*(randn(M, T) + 1i * randn(M, T));
    Ek = Nk'*A_R;
    Yk = Phi * Xk + Ek;
end
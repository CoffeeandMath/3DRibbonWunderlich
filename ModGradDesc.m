function q = ModGradDesc(f,df,ddf,q0)
tol = 10^-8;

fnorm = 10;
kmax = 100;
q = q0;

k = 0;
while and(k<kmax,fnorm > tol)
    k = k+1;
    if mod(k,10)==1
        Amat = ddf(q);
    end
    
    q = q - Amat\df(q);
    
    fnorm = df(q);
end
end


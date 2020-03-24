function [Xmin] = nlinconjgrad(f,Df,x0,alpha)
N = max(size(x0));
tol = 10^-7;
dfk = Df(x0);
Nalpha = 50;
p = - dfk;
k = 0;
alphai = linspace(alpha/Nalpha,alpha,Nalpha);
xk = x0;

while(norm(dfk)>tol)
    %alpha = graddesc(@(gamma) f(xk+gamma*p),@(gamma) Df(xk+gamma*p),0);
    fnorm = norm(dfk);
    fmin = zeros(Nalpha,1);
    for i = 1:Nalpha
        xkpi = xk + alphai(i)*p;
        fmin(i) = f(xkpi);
        if fmin(i) > 10^6 || isnan(fmin(i)) || imag(fmin(i)) ~= 0
            fmin(i) = Inf;
        end
    end
    [~,I] = min(fmin);
    alphaimin = alphai(I);
    if I == 1
        alphain = linspace(alpha/(Nalpha^2),alpha/Nalpha,Nalpha);
        p = - Df(xk);
        for i = 1:Nalpha
            xkpi = xk + alphain(i)*p;
            fmin(i) = f(xkpi);
            if fmin(i) > 10^6 || isnan(fmin(i)) || imag(fmin(i)) > 0.01
                fmin(i) = inf;
            end
        end
        [~,I] = min(fmin);
        alphaimin = alphain(I);
    end
    
    xkp  = xk + alphaimin*p;
    dfkp = Df(xkp);
    
    beta = norm(dfkp)^2/norm(dfk)^2;
    p = - dfkp + beta*p;
    k = k+1;
    xk = xkp;
    dfk = dfkp;
end
Xmin = xk;

k
end


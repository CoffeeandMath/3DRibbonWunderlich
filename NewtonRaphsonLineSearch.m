function [qout,k] = NewtonRaphsonLineSearch(f,df,A,q0)

%tol = 10^-8 * dT*sqrt(D/m)/dS;
%dfnorm = norm(df(q0));
dfnorm = 10;
dim = max(size(q0));
%tol = min([fnorm/10^3,10^-8 * dT*sqrt(D/mtot)/dS]);
tol = sqrt(dim)*10^-6;
kmax = 100;
qn = q0;
k = 0;
% Nls = 10;
% alphai = linspace(-.1,1.5,Nls);
Nls = 1;
alphai = 1;
while and(dfnorm > tol,k < kmax)
    k = k+1;
    Amat = A(qn);
    
    stepdir = -( (Amat\df(qn)));
    if or(norm(1.0*isnan(stepdir))>0,norm(stepdir)==Inf)
        stepdir = -df(qn);
    end
    fval = zeros(Nls,1);
    for i = 1:Nls
        fval(i) = f(qn+alphai(i)*stepdir);
        
        if or(fval(i) > 10^8,isnan(fval(i)))
            fval(i) = Inf;
        end
    end
    [~,I] = min(fval);
%     if I == Nls
%         fprintf('railed')
%     end
    step = alphai(I)*stepdir;
    
    
    
    stepn = norm(step);
    while norm(1.0*isnan(df(qn+step))) > 1
        step = step/2;
    end
    qn = qn + step;
    dfnorm = norm(df(qn));
end
%k
qout = qn;
end
function [qout,k] = NewtonRaphson(f,df,A,q0)
global Nn;
%tol = 10^-8 * dT*sqrt(D/m)/dS;
dfnorm = norm(df(q0));
%tol = min([fnorm/10^3,10^-8 * dT*sqrt(D/mtot)/dS]);
tol = Nn*10^-12;
kmax = 100;
qn = q0;
k = 0;
Ndof = max(size(q0));

while and(dfnorm > tol,k < kmax)
    k = k+1;
    stepdir = -( (A(qn)\df(qn)));
   
    qn = qn + stepdir;
    dfnorm = norm(df(qn));
end
k
qout = qn;
end
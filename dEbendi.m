function dEout = dEbendi(kappa1,kappa3,Sys)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Nn = Sys.Nn;
dS = Sys.dS;
[omega,eta] = kappa2omegaeta(kappa1,kappa3,dS);

D = Sys.D;
dEout = zeros(Nn,3);

dEout(:,1) = D(:,1).*omega.*(1-eta.^4);
dEout(:,3) = D(:,1).*2.*omega.*eta.*(1+eta.^2);

end


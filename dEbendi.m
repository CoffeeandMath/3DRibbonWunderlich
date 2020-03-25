function dEout = dEbendi(thetai,Sys)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Nn = Sys.Nn;
dS = Sys.dS;
[omega,eta] = kappa2omegaeta(thetai(:,1),thetai(:,3),dS);

D = Sys.D;
dEout = zeros(Nn,3);

dEout(:,1) = D(:,1).*omega.*(1-eta.^4);
dEout(:,2) = D(:,2).*thetai(:,2);
dEout(:,3) = D(:,1).*2.*omega.*eta.*(1+eta.^2);

end


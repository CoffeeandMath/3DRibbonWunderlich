function Eout = Ebendi(thetai,Sys)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dS = Sys.dS;
D = Sys.D;
[omega,eta] = kappa2omegaeta(thetai(:,1),thetai(:,3),dS);
Eout = 0.5*D(:,1).*omega.^2.*(1+eta.^2).^2 + 0.5*D(:,2).*(thetai(:,2).^2);

end


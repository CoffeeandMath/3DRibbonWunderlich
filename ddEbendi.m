function ddEbendiout = ddEbendi(thetai,Sys)
Nn = Sys.Nn;
D = Sys.D;
dS = Sys.dS;
[~,eta] = kappa2omegaeta(thetai(:,1),thetai(:,3),dS);
ddEbendiout = zeros(3,3,Nn);
for i = 2:(Nn-1)

    et = eta(i);
    ddEbendiout(:,:,i) = [D(i,1)*(1+3*et^4), 0, -D(i,1)*et^3;...
        0, D(i,2), 0;...
        -D(i,1)*et^3, 0, 2*D(i,1)*(1 + 3*et^2)];
end

end


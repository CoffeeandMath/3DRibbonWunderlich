function dkout = ddkappadphidphi(qi)
dkout = zeros(3,2,2);

dkout(:,1,1) = -(1/2)*[0,1,0,0;0,0,1,0;0,0,0,1]*qi;
dkout(:,1,2) = -(1/2)*[0,1,0,0;0,0,1,0;0,0,0,1]*qmult([0;0;0;1],qi,[0;0;0;1]);
dkout(:,2,1) = dkout(:,1,2);
dkout(:,2,2) = dkout(:,1,1);
end


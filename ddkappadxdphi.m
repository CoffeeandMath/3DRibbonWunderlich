function ddkout = ddkappadxdphi(philoc,tloc,dloc,Tloc,Dloc,ptwloc,pparloc,epsvalloc)
%% Setting local values
phinm = philoc{1};
phin = philoc{2};

tnm = tloc{1};
tn = tloc{2};

ptwnm = ptwloc{1};
ptwn = ptwloc{2};

pparnm = pparloc{1};
pparn = pparloc{2};

Dnm = Dloc{1};
Dn = Dloc{2};

dnm = dloc{1};
dn = dloc{2};

Tnm =  Tloc{1};
Tn = Tloc{2};

epsvalnm = epsvalloc{1};
epsvaln = epsvalloc{2};
%%





I3 = [0,1,0,0;0,0,1,0;0,0,0,1];
skw = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
e3 = [0;0;0;1];

fnm = sqrt(1+Tnm'*tnm);
Arnm = (0.5/fnm)*((-1/fnm)*(cross(Tnm,tnm)*Tnm') + fnm * skw(Tnm) - (1/fnm)*(Tnm*cross(tnm,Tnm)'));
Brnm = Arnm*(1/epsvalnm)*(eye(3) - tnm*tnm');
dqndphinmdxnm = -[qmult(e3,qconj(dnm),[0;Brnm(:,1)],dn),qmult(e3,qconj(dnm),[0;Brnm(:,2)],dn),qmult(e3,qconj(dnm),[0;Brnm(:,3)],dn)];
ddkout(:,:,1,1) = I3*dqndphinmdxnm;


fn = sqrt(1+Tn'*tn);
Arn = (0.5/fn)*((-1/fn)*(cross(Tn,tn)*Tn') + fn * skw(Tn) - (1/fn)*(Tn*cross(tn,Tn)'));
Brn = Arn*(1/epsvaln)*(eye(3) - tn*tn');
dqndphinmdxnp = -[qmult(e3,qconj(dnm),[0;Brn(:,1)],dn),qmult(e3,qconj(dnm),[0;Brn(:,2)],dn),qmult(e3,qconj(dnm),[0;Brn(:,3)],dn)];
ddkout(:,:,3,1) = I3*dqndphinmdxnp;

ddkout(:,:,2,1) = -ddkout(:,:,1,1) - ddkout(:,:,3,1);



dqndphindxnm = [qmult(qconj(dnm),[0;Brnm(:,1)],dn,e3),qmult(qconj(dnm),[0;Brnm(:,2)],dn,e3),qmult(qconj(dnm),[0;Brnm(:,3)],dn,e3)];
ddkout(:,:,1,2) = I3*dqndphindxnm;

dqndphindxnp = [qmult(qconj(dnm),[0;Brn(:,1)],dn,e3),qmult(qconj(dnm),[0;Brn(:,2)],dn,e3),qmult(qconj(dnm),[0;Brn(:,3)],dn,e3)];
ddkout(:,:,3,2) = I3*dqndphindxnp;

ddkout(:,:,2,2) = -ddkout(:,:,1,2) - ddkout(:,:,3,2);

end


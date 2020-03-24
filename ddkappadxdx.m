function ddkout = ddkappadxdx(tloc,dloc,Tloc,Dloc,ptwloc,epsvalloc)


%% ddkappaIdxIndxKm
%upper case are components, lower case are nodes
% 1)component of kappa
% 2)component of x
% 3)component of x
% 4)node x
% 5)node x
ddkout = zeros(3,3,3,3,3);

skw = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];


tim = tloc{1};
ti = tloc{2};
Tim = Tloc{1};
Ti = Tloc{2};
epsim = epsvalloc{1};
epsi = epsvalloc{2};
dim = dloc{1};
di = dloc{2};

f2im = 1 + Tim'*tim;
fim = sqrt(f2im);

f2i = 1 + Ti'*ti;
fi = sqrt(f2i);

Aim = (0.5/fim)*((-1/fim)*(cross(Tim,tim)*Tim') + fim * skw(Tim) - (1/fim)*(Tim*cross(tim,Tim)'));
Bim = Aim*(1/epsim)*(eye(3) - tim*tim');


Ai = (0.5/fi)*((-1/fi)*(cross(Ti,ti)*Ti') + fi * skw(Ti) - (1/fi)*(Ti*cross(ti,Ti)'));
Bi = Ai*(1/epsi)*(eye(3) - ti*ti');


dBimdeps = thirdorder(Bim,-tim/epsim);
Timctim = cross(Tim,tim);
dBimdf2 = thirdorder(Tim*Timctim' - Timctim*Tim' - Timctim*tim',(0.5/epsim^2)*(-1/f2im^2)*(Tim-tim*(f2im-1)));
Timskw = skw(Tim);
timprojeps = eye(3) - tim*tim';
dBimdtim = (0.5/(epsim^2*f2im))*(thirdorder(Tim,Timskw - Timctim*tim') - thirdorderin(Timskw - Timctim*tim',tim+Tim) - thirdorder(Timctim,timprojeps));
dBim = dBimdeps + dBimdf2 + dBimdtim;

dBideps = thirdorder(Bi,-ti/epsi);
Ticti = cross(Ti,ti);
dBidf2 = thirdorder(Ti*Ticti' - Ticti*Ti' - Ticti*ti',(0.5/epsi^2)*(-1/f2i^2)*(Ti-ti*(f2i-1)));
Tiskw = skw(Ti);
tiprojeps = eye(3) - ti*ti';
dBidti = (0.5/(epsi^2*f2i))*(thirdorder(Ti,Tiskw - Ticti*ti') - thirdorderin(Tiskw - Ticti*ti',ti+Ti) - thirdorder(Ticti,tiprojeps));
dBi = dBideps + dBidf2 + dBidti;


dBimBimd = zeros(4,3,3);
dBimBid = zeros(4,3,3);
dBiBid = zeros(4,3,3);

cdim = qconj(dim);

for n = 1:3
    for q = 1:3
        
        dBimBimd(:,n,q) = qmult(cdim,qmult([0;Bim(:,q)],[0;Bim(:,n)]) - [0;dBim(:,n,q)],di);
        dBimBid(:,n,q) = qmult(cdim,[0;Bim(:,n)],[0;Bi(:,q)],di);
        dBiBid(:,n,q) = qmult(cdim,qmult([0;Bi(:,n)],[0;Bi(:,q)]) + [0;dBi(:,n,q)],di);
    end
end
ddkout(:,:,:,1,1) = 2*dBimBimd(2:4,:,:);
ddkout(:,:,:,1,3) = 2*dBimBid(2:4,:,:);
ddkout(:,:,:,3,3) = 2*dBiBid(2:4,:,:);
ddkout(:,:,:,3,1) = permute(ddkout(:,:,:,1,3),[1,3,2]);
ddkout(:,:,:,1,2) = -ddkout(:,:,:,1,1) - ddkout(:,:,:,1,3);
ddkout(:,:,:,2,1) = permute(ddkout(:,:,:,1,2),[1,3,2]);
ddkout(:,:,:,2,3) = -ddkout(:,:,:,1,3) - ddkout(:,:,:,3,3);
ddkout(:,:,:,3,2) = permute(ddkout(:,:,:,2,3),[1,3,2]);
ddkout(:,:,:,2,2) = -ddkout(:,:,:,2,3) - ddkout(:,:,:,2,1);



end


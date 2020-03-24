 function dkappadxout = dkappadx(epsim,epsi,tim,ti,Tim,Ti,Dim,Di,dim,di,ptwim,ptwi)
dkappadxout = zeros(3,3,3);

skw = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];

fim = sqrt(1+Tim'*tim);
Arim = (0.5/fim)*((-1/fim)*(cross(Tim,tim)*Tim') + fim * skw(Tim) - (1/fim)*(Tim*cross(tim,Tim)'));
Brim = Arim*(1/epsim)*(eye(3) - tim*tim');
dqidxim = [qmult(qconj(dim),[0;Brim(:,1)],di),qmult(qconj(dim),[0;Brim(:,2)],di),qmult(qconj(dim),[0;Brim(:,3)],di)];
dkappadxout(:,:,1) = 2*dqidxim(2:4,:);


fi = sqrt(1+Ti'*ti);
Ari = (0.5/fi)*((-1/fi)*(cross(Ti,ti)*Ti') + fi * skw(Ti) - (1/fi)*(Ti*cross(ti,Ti)'));
Bri = Ari*(1/epsi)*(eye(3) - ti*ti');
dqidxip = [qmult(qconj(dim),[0;Bri(:,1)],di),qmult(qconj(dim),[0;Bri(:,2)],di),qmult(qconj(dim),[0;Bri(:,3)],di)];
dkappadxout(:,:,3) = 2*dqidxip(2:4,:);

dkappadxout(:,:,2) = -dkappadxout(:,:,1) - dkappadxout(:,:,3);


end


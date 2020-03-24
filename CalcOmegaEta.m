function [omega,eta] = CalcOmegaEta(X,Sys)


Nn = Sys.Nn;
dS = Sys.dS;
dstar = Sys.dstar;

x = cell(Nn,1);
phi = zeros(Nn-1,1);

for i = 1:Nn
    x{i} = X((4*i - 3):(4*i - 1));
    if i < Nn
        phi(i) = X(4*i);
    end
end



%% Calculating the tangents
dcur = cell(Nn-1,3);

for i = 1:(Nn-1)
    dcur{i,3} = (x{i+1} - x{i})/norm(x{i+1} - x{i});
end


%% Calculating rotation from cartesian to reference frames
Di = cell(Nn,1);
% for i = 1:(Nn-1)
%     pdIdoteI = 1 + dstar{i,1}'*[1;0;0] + dstar{i,2}'*[0;1;0] + dstar{i,3}'*[0;0;1];
%     epsR = cross(dstar{i,1},[1;0;0]) + cross(dstar{i,2},[0;1;0]) + cross(dstar{i,3},[0;0;1]);
%     r = (1/2) * sqrt(pdIdoteI);
%     v = -(1/2) * epsR/sqrt(pdIdoteI);
%     Di{i} = qgen(r,v(1),v(2),v(3));
% end
i = 1;
pdIdoteI = 1 + dstar{i,1}'*[1;0;0] + dstar{i,2}'*[0;1;0] + dstar{i,3}'*[0;0;1];
epsR = cross(dstar{i,1},[1;0;0]) + cross(dstar{i,2},[0;1;0]) + cross(dstar{i,3},[0;0;1]);
r = (1/2) * sqrt(pdIdoteI);
v = -(1/2) * epsR/sqrt(pdIdoteI);
Di{i} = qgen(r,v(1),v(2),v(3));

for i = 2:(Nn-1)
    pdIdoteI = 1 + dstar{i,1}'*dstar{i-1,1} + dstar{i,2}'*dstar{i-1,2} + dstar{i,3}'*dstar{i-1,3};
    epsR = cross(dstar{i,1},dstar{i-1,1}) + cross(dstar{i,2},dstar{i-1,2}) + cross(dstar{i,3},dstar{i-1,3});
    r = (1/2) * sqrt(pdIdoteI);
    v = -(1/2) * epsR/sqrt(pdIdoteI);
    Direl = qgen(r,v);
    Di{i} = qmult(Direl,Di{i-1});
end
%% Calculating the twist quaternion
ptwist = cell(Nn-1,1);

for i = 1:(Nn-1)
    r = cos(phi(i)/2);
    v = sin(phi(i)/2)*dstar{i,3};
    ptwist{i} = qgen(r,v(1),v(2),v(3));
end

%% Calculating parallel transport quaternion
ppar = cell(Nn-1,1);

for i = 1:(Nn-1)
    
    r = sqrt((1+dstar{i,3}'*dcur{i,3})/2);
    v = (1/sqrt(2))*(1/sqrt(1 + dstar{i,3}'*dcur{i,3})) * cross(dstar{i,3},dcur{i,3});
    ppar{i} = qgen(r,v(1),v(2),v(3));
end
%% Calculating the product of quaternions
di = cell(Nn-1,1);

for i = 1:(Nn-1)
    di{i} = qmult(ppar{i},ptwist{i},Di{i});
end


%% Calculating Relative Rotation quaternion
q_i = cell(Nn,1);

for i = 2:(Nn-1)
    q_i{i} = qmult(qconj(di{i-1}),di{i});
end

%% Calculating the axial quaternion

kappa = cell(Nn,1);
for i = 2:(Nn-1)
    %kappa{i} = 2*qlog(q_i{i});
    kappa{i} = 2*q_i{i};
%     [~,pv1,pv2,pv3] = parts(p_i{i});
%     theta{i} = quaternion(0,pv1,pv2,pv3);
end


%% Calculating strain measures
thetai = zeros(Nn,3);

for i = 2:(Nn-1)

    thetai(i,1) = kappa{i}(2);
    thetai(i,2) = kappa{i}(3);
    thetai(i,3) = kappa{i}(4);
end

[omega,eta] = kappa2omegaeta(thetai(:,1),thetai(:,3),dS);
end


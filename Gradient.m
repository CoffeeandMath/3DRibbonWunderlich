function dE = Gradient(X,Sys)


Nn = Sys.Nn;
D = Sys.D;
dS = Sys.dS;
rho = Sys.rho;
Eaxial = Sys.Eaxial;
Fxnode = Sys.Fxnode;
Fynode = Sys.Fynode;
Fznode = Sys.Fznode;
dstar = Sys.dstar;

x = cell(Nn,1);
for i = 1:Nn
    x{i} = X((4*i - 3):(4*i - 1));
end
phi = X(4:4:end);



%% Calculating the tangents
dcur = cell(Nn-1,3);
epsval = zeros(Nn-1,1);
for i = 1:(Nn-1)
    epsval(i) = norm(x{i+1} - x{i});
    dcur{i,3} = (x{i+1} - x{i})/epsval(i);
end


%% Calculating rotation from cartesian to reference frames
Di = cell(Nn,1);

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
    ptwist{i} = qgen(r,v);
end

%% Calculating parallel transport quaternion
ppar = cell(Nn-1,1);

for i = 1:(Nn-1)
    
    r = sqrt((1+dstar{i,3}'*dcur{i,3})/2);
    v = (1/sqrt(2))*(1/sqrt(1 + dstar{i,3}'*dcur{i,3})) * cross(dstar{i,3},dcur{i,3});
    ppar{i} = qgen(r,v);
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

%% Now to Calculate The Gradients
Bi = zeros(3,3,Nn-1);
skw = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
for i =1:(Nn-1)
    Ti = dstar{i,3};
    ti = dcur{i,3};
    Bi(:,:,i) = (0.5/epsval(i))*((1/(1+Ti'*ti))*(Ti*cross(Ti,ti)' - cross(Ti,ti)*((ti+Ti)')) + skw(Ti));
end


dkappadphivals = zeros(3,2,Nn-1);
dkappadxvals = zeros(3,3,3,Nn); 

for i = 2:(Nn-1)
    dkappadxvals(:,:,:,i) = dkappadx(di{i-1},di{i},Bi(:,:,i-1),Bi(:,:,i));
    dkappadphivals(:,:,i) = dkappadphi(q_i{i}); 
    
end


%% Calculating derivatives of energy

dEbenddgamma = zeros(Nn-1,1);
dEbenddkappa = dEbendi(thetai,Sys);
for i = 2:(Nn-1)
    
    dEbenddgamma((i-1):i) = dEbenddgamma((i-1):i) + (dEbenddkappa(i,:)*dkappadphivals(:,:,i))';
    
end

dEbenddx1 = zeros(Nn,1);
dEbenddx2 = zeros(Nn,1);
dEbenddx3 = zeros(Nn,1);

for i = 2:(Nn-1)

%     dEbenddx1((i-1):(i+1)) = dEbenddx1((i-1):(i+1)) + (dEbenddkappa(i,1)*dkappa1dx{1}(i,:) + dEbenddkappa(i,2)*dkappa2dx{1}(i,:) + dEbenddkappa(i,3)*dkappa3dx{1}(i,:))';
%     dEbenddx2((i-1):(i+1)) = dEbenddx2((i-1):(i+1)) + (dEbenddkappa(i,1)*dkappa1dx{2}(i,:) + dEbenddkappa(i,2)*dkappa2dx{2}(i,:) + dEbenddkappa(i,3)*dkappa3dx{2}(i,:))';
%     dEbenddx3((i-1):(i+1)) = dEbenddx3((i-1):(i+1)) + (dEbenddkappa(i,1)*dkappa1dx{3}(i,:) + dEbenddkappa(i,2)*dkappa2dx{3}(i,:) + dEbenddkappa(i,3)*dkappa3dx{3}(i,:))';
    
    dEbenddx1((i-1):(i+1)) = dEbenddx1((i-1):(i+1)) + (dEbenddkappa(i,:)*squeeze(dkappadxvals(:,1,:,i)))';
    dEbenddx2((i-1):(i+1)) = dEbenddx2((i-1):(i+1)) + (dEbenddkappa(i,:)*squeeze(dkappadxvals(:,2,:,i)))';
    dEbenddx3((i-1):(i+1)) = dEbenddx3((i-1):(i+1)) + (dEbenddkappa(i,:)*squeeze(dkappadxvals(:,3,:,i)))';
end


dEgravdgamma = zeros(Nn-1,1);

dEgravdx1 = zeros(Nn,1);
dEgravdx2 = zeros(Nn,1);
dEgravdx3 = rho*ones(Nn,1);


dEaxdgamma = zeros(Nn-1,1);
dEaxdx1 = zeros(Nn,1);
dEaxdx2 = zeros(Nn,1);
dEaxdx3 = zeros(Nn,1);


for i = 1:(Nn)
    if i == 1
      
        v = -Eaxial*(epsval(i) - dS)*dcur{i,3};
        dEaxdx1(i) = v(1);
        dEaxdx2(i) = v(2);
        dEaxdx3(i) = v(3);
    elseif i == Nn
       
        v = Eaxial*(epsval(i-1) - dS)*dcur{i-1,3};
        dEaxdx1(i) = v(1);
        dEaxdx2(i) = v(2);
        dEaxdx3(i) = v(3);
    else
        v = Eaxial*(epsval(i-1) - dS)*dcur{i-1,3} - Eaxial*(epsval(i) - dS)*dcur{i,3};
        dEaxdx1(i) = v(1);
        dEaxdx2(i) = v(2);
        dEaxdx3(i) = v(3);
    end
end


dEdgamma = dEbenddgamma + dEgravdgamma + dEaxdgamma;
dEdx1 = dEbenddx1 + dEgravdx1 + dEaxdx1 - Fxnode;
dEdx2 = dEbenddx2 + dEgravdx2 + dEaxdx2 - Fynode;
dEdx3 = dEbenddx3 + dEgravdx3 + dEaxdx3 - Fznode;

dE = zeros(4*Nn-1,1);
for i = 1:Nn
    dE(4*i-3) = dEdx1(i);
    dE(4*i-2) = dEdx2(i);
    dE(4*i-1) = dEdx3(i);
    if i < Nn
        dE(4*i) = dEdgamma(i);
    end
end


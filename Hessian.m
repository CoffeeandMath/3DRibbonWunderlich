function H = Hessian(X,Sys)


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
phi = zeros(Nn-1,1);

for i = 1:Nn
    x{i} = X((4*i - 3):(4*i - 1));
    if i < Nn
        phi(i) = X(4*i);
    end
end



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

dkappadphivals = zeros(3,2,Nn-1);



dkappadxvals = zeros(3,3,3,Nn);

for i = 2:(Nn-1)

    dkappadxvals(:,:,:,i) = dkappadx(epsval(i-1),epsval(i),dcur{i-1,3},dcur{i,3},dstar{i-1,3},dstar{i,3},Di{i-1},Di{i},di{i-1},di{i},ptwist{i-1},ptwist{i});
    dkappadphivals(:,:,i) = dkappadphi(q_i{i});

    
end

%% Grouping parameters to avoid parallel broadcast variable issues
philoc = cell(Nn-1,1);
tloc = cell(Nn-1,1);
diloc = cell(Nn-1,1);
Tloc = cell(Nn-1,1);
Diloc = cell(Nn-1,1);
ptwloc = cell(Nn-1,1);
pparloc = cell(Nn-1,1);
epsvalloc = cell(Nn-1,1);

%(phinm,phin,tnm,tn,dnm,dn,Tnm,Tn,Dnm,Dn,ptwnm,ptwn,pparnm,pparn,epsvalnm,epsvaln)
for i = 2:(Nn-1)
    philoc{i} = cell(2,1);
    philoc{i}{1} = phi(i-1);
    philoc{i}{2} = phi(i);
    
    tloc{i} = cell(2,1);
    tloc{i}{1} = dcur{i-1,3};
    tloc{i}{2} = dcur{i,3};
    
    diloc{i} = cell(2,1);
    diloc{i}{1} = di{i-1};
    diloc{i}{2} = di{i};
    
    Tloc{i} = cell(2,1);
    Tloc{i}{1} = dstar{i-1,3};
    Tloc{i}{2} =  dstar{i,3};
    
    Diloc{i} = cell(2,1);
    Diloc{i}{1} = Di{i-1};
    Diloc{i}{2} = Di{i};
    
    ptwloc{i} = cell(2,1);
    ptwloc{i}{1} = ptwist{i-1};
    ptwloc{i}{2} = ptwist{i};
    
    pparloc{i} = cell(2,1);
    pparloc{i}{1} = ppar{i-1};
    pparloc{i}{2} = ppar{i};
    
    epsvalloc{i} = cell(2,1);
    epsvalloc{i}{1} = epsval(i-1);
    epsvalloc{i}{2} = epsval(i);
    
end

%% Calculating second derivatives of strain measures



%ddkappaIdxkJdphiK

ddkappadphidphival = zeros(3,2,2,Nn); %ddkappaiNdphijdphiK = a(i,j,K,N)
ddkappadxdphival = zeros(3,3,3,2,Nn); %ddkappaiNdxjPdphiQ = a(i,j,P,Q,N)
ddkappadxdxval = zeros(3,3,3,3,3,Nn); %ddkappaiNdxjPdxkQ = a(i,j,k,P,Q,N)

ddEbend = spalloc(4*Nn-1,4*Nn-1,20*Nn);
%ddEbend = zeros(4*Nn-1,4*Nn-1);
ddEbendival = ddEbendi(thetai,Sys);
dEbendival = dEbendi(thetai,Sys);
ddEbendlocal = cell(Nn,1);


parfor n = 2:(Nn-1)
    
    %dp hinm dphinm
    
    
    ddkappadphidphival(:,:,:,n) = ddkappadphidphi(q_i{n});
    

    %ddkappaIdxkJdphiK
    ddkappadxdphival(:,:,:,:,n) = ddkappadxdphi(philoc{n},tloc{n},diloc{n},Tloc{n},Diloc{n},ptwloc{n},pparloc{n},epsvalloc{n});

    %ddkappaIdxkJdxmL
    
    ddkappadxdxval(:,:,:,:,:,n) = ddkappadxdx(tloc{n},diloc{n},Tloc{n},Diloc{n},ptwloc{n},epsvalloc{n});
end



llocal = zeros(11,11);
for j = 2:(Nn-1)
    
    
    
    %dximdxim
    llocal(1:3,1:3) = dkappadxvals(:,:,1,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,1,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,1,1,j));
    %dximdxi
    llocal(1:3,5:7) = dkappadxvals(:,:,1,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,2,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,1,2,j));
    %dxidxim
    llocal(5:7,1:3) = llocal(1:3,5:7)';
    %dximdxip
    llocal(1:3,9:11) = dkappadxvals(:,:,1,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,3,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,1,3,j));
    %dxipdxim
    llocal(9:11,1:3) = llocal(1:3,9:11)';
    %dxidxi
    llocal(5:7,5:7) = dkappadxvals(:,:,2,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,2,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,2,2,j));
    %dxidxip
    llocal(5:7,9:11) = dkappadxvals(:,:,2,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,3,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,2,3,j));
    %dxipdxi
    llocal(9:11,5:7) = llocal(5:7,9:11)';
    %dxipdxip
    llocal(9:11,9:11) = dkappadxvals(:,:,3,j)'*ddEbendival(:,:,j)*dkappadxvals(:,:,3,j) + vtimes3D(dEbendival(j,:),ddkappadxdxval(:,:,:,3,3,j));
    
    
    
    %dximdphi
    llocal(1:3,[4,8]) = dkappadxvals(:,:,1,j)'*ddEbendival(:,:,j)*dkappadphivals(:,:,j) + vtimes3D(dEbendival(j,:),squeeze(ddkappadxdphival(:,:,1,:,j)));
    llocal([4,8],1:3) = llocal(1:3,[4,8])';
    %dxidphi
    llocal(5:7,[4,8]) = dkappadxvals(:,:,2,j)'*ddEbendival(:,:,j)*dkappadphivals(:,:,j) + vtimes3D(dEbendival(j,:),squeeze(ddkappadxdphival(:,:,2,:,j)));
    llocal([4,8],5:7) = llocal(5:7,[4,8])';
    %dxipdphi
    llocal(9:11,[4,8]) = dkappadxvals(:,:,3,j)'*ddEbendival(:,:,j)*dkappadphivals(:,:,j) + vtimes3D(dEbendival(j,:),squeeze(ddkappadxdphival(:,:,3,:,j)));
    llocal([4,8],9:11) = llocal(9:11,[4,8])';
    
    
    llocal([4,8],[4,8]) = dkappadphivals(:,:,j)'*ddEbendival(:,:,j)*dkappadphivals(:,:,j) + vtimes3D(dEbendival(j,:),ddkappadphidphival(:,:,:,j));
    
    
    ddEbendlocal{j} = llocal;
end



for i = 2:(Nn-1)
    ddEbend(((4*i-7):(4*i+3)),((4*i-7):(4*i+3))) = ddEbend(((4*i-7):(4*i+3)),((4*i-7):(4*i+3))) + ddEbendlocal{i};
end


%% Eaxial hessian
ddEaxloc = cell(Nn-1,1);
ddEax = spalloc(4*Nn-1,4*Nn-1,9*Nn);
for i = 1:(Nn-1)
    ddEaxloc{i} = ddestrdx(dcur{i,3},epsval(i),dS);
end

for i = 1:(Nn-1)
    ddEax((4*i-3):(4*i-1),(4*i-3):(4*i-1)) = ddEax((4*i-3):(4*i-1),(4*i-3):(4*i-1)) + Eaxial*ddEaxloc{i}(1:3,1:3);
    ddEax((4*i-3):(4*i-1),(4*i+1):(4*i+3)) = ddEax((4*i-3):(4*i-1),(4*i+1):(4*i+3)) + Eaxial*ddEaxloc{i}(1:3,4:6);
    ddEax((4*i+1):(4*i+3),(4*i-3):(4*i-1)) = ddEax((4*i+1):(4*i+3),(4*i-3):(4*i-1)) + Eaxial*ddEaxloc{i}(4:6,1:3);
    ddEax((4*i+1):(4*i+3),(4*i+1):(4*i+3)) = ddEax((4*i+1):(4*i+3),(4*i+1):(4*i+3)) + Eaxial*ddEaxloc{i}(4:6,4:6);
end




%% Putting it together
H = ddEbend + ddEax;

end
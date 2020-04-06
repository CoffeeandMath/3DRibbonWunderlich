function Etot = Objective(X,Sys)
% We want to write a script which returns the objective for optimization.
% this function takes in the input variables X which are a columnated
% x1,x2,x3,gamma where xi are the node positions and gamma is a parallel
% transport twist for the frames. dstar is a cell array (size Nn-1,3) at the previous
% time step. Note that X has dim 4*Nn-1 because there are only Nn-1
% segments and gamma is defined at the segments
Nn = Sys.Nn;
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

for i = 1:(Nn-1)
    dcur{i,3} = (x{i+1} - x{i})/norm(x{i+1} - x{i});
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
%% Calculating strain energy
% Ebend = 0;
% 
% for i = 2:(Nn-1)
%     Ebend = Ebend + (1/2)*(dot(D(i,:),thetai(i,:).^2));
% end
Ebend = Ebendi(thetai,Sys);

%% Calculating gravitational potential energy


Egrav = 0;

for i = 1:Nn
    Egrav = Egrav + rho*x{i}(3);
end

%% Calculating axial strain
Eax = 0;
for i = 1:(Nn-1)
    Eax = Eax + (1/2)*Eaxial*(norm(x{i+1} - x{i}) - dS)^2;
end

%% Potential energy of externally applied forces
Eext = 0;
for i = 1:Nn
    Eext = Eext + Fxnode(i)*x{i}(1) + Fynode(i)*x{i}(2) + Fznode(i)*x{i}(3);
end
%% Calculating Total Energy
Etot = Ebend + Egrav + Eax - Eext;


%% Auxilatory Functions


end


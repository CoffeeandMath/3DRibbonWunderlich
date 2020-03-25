%%
clear all
close all
ploton = false;
savedata = true;
%%


Sys = struct();

Nn = 120;
L = 1;
w = L/10;
Si = linspace(0,L,Nn);
dS = Si(2) - Si(1);

rho = 0*10^-1/Nn;
Eaxial = 1000/dS;
D = 2*10^-2*ones(Nn,2)/dS; D(:,2) = 1000*D(:,2);


Fxnode = zeros(Nn,1); Fxnode(end) = -0.0001;
Fynode = zeros(Nn,1); Fynode(end) = .5;
Fznode = zeros(Nn,1); Fznode(end) = .3;



sc = L/30;

dstar0 = cell(Nn-1,3);
d10 = [0;1;0];
d20 = [0;0;1];
d30 = [1;0;0];

for i = 1:(Nn-1)

    thsc = -(i-1)/(Nn*100);
    R = [cos(thsc), 0 , sin(thsc);...
        0, 1, 0;...
        -sin(thsc), 0, cos(thsc)];
    dstar0{i,1} = R*d10;
    dstar0{i,2} = R*d20;
    dstar0{i,3} = R*d30;
end

xinit = cell(Nn,1);
xinit{1} = zeros(3,1);
for i = 1:(Nn-1)
    xinit{i+1} = xinit{i} + dS*dstar0{i,3};
    
end

gammainit = zeros(Nn-1,1);

Xinit = zeros(4*Nn-1,1);
for i = 1:Nn
    Xinit((4*i-3):(4*i-1)) = xinit{i} ;
    if i < Nn
        Xinit(4*i) = gammainit(i);
    end
end
Xinitconst = Xinit(8:end);



rho0 = rho;
Nk = 1500;
mult = linspace(0,1,Nk);

dstar = dstar0;
figure


Sys.Nn = Nn;
Sys.D = D;
Sys.rho = rho;
Sys.dS = dS;
Sys.Eaxial = Eaxial;
Sys.Fxnode = Fxnode;
Sys.Fynode = Fynode;
Sys.Fznode = Fznode;
Sys.dstar = dstar;
for i = 1:Nk
    if mod(i,10)==0
        fprintf(['step = ' num2str(i) ', ' num2str(i*100/Nk) ' percent done \n'])
    end
    %tic;
    rho = mult(i)*rho0;
    Sys.rho = rho;
    Sys.Fxnode = mult(i)*Fxnode;
    Sys.Fynode = mult(i)*Fynode;
    Sys.Fznode = mult(i)*Fznode;
    Sys.dstar = dstar;
    ff.objective = @(X) Objective(Xextender(X,dS),Sys);
    %ff.gradient = @(X) CD(ff.objective,X);
    ff.gradient = @(X) Xconstrainer(Gradient(Xextender(X,dS),Sys));
    %ff.gradient2 = @(X) Xconstrainer(Gradient2(Xextender(X),dstar));
    %ff.hessian = @(X) CDGradient(ff.gradient,X);
    ff.hessian = @(X) Hconstrainer(Hessian(Xextender(X,dS),Sys));
    %ff.hessian = @(X) Hconstrainer(Hessian(Xextender(X,dS),Sys));
    options = optimoptions('fminunc');
    options.MaxIterations = 1000;
    options.StepTolerance = 10^-10;
    options.OptimalityTolerance = 10^-10;
    %Xsolve = nlinconjgrad(ff.objective,ff.gradient,Xinitconst,10^-7);
    %Xsolve = fminunc(ff.objective,Xinitconst,options);
    %Xsolve = NewtonRaphsonLineSearch(ff.objective,ff.gradient,ff.hessian,Xinitconst);
    %norm(ff.gradient2(Xsolve))
    Xsolve = NewtonRaphsonLineSearch(ff.objective,ff.gradient,ff.hessian,Xinitconst);
    
    x = cell(Nn,1);
    x{1} = [0;0;0];
    x{2} = dS*[1;0;0];
    gamma = zeros(Nn-1,1);
    for j = 3:(Nn)
        x{j} = Xsolve((4*j-10):(4*j-8));
    end
    for j = 2:(Nn-1)
        gamma(j) = Xsolve(4*j-7);
    end
    Xinitconst = Xsolve;
    %Xinitconst(1:4:end) = 0;
    Xvals = zeros(Nn,3);
    Xseg = zeros(Nn-1,3);
    for j = 1:Nn
        Xvals(j,:) = x{j}';
        if j < Nn
            Xseg(j,:) = (x{j}+ x{j+1})/2';
        end
    end
    dcur = cell(Nn-1,3);
    for j = 1:(Nn-1)
        
        dcur{j,3} = (x{j+1} - x{j})/norm(x{j+1} - x{j});
    end
    %% Calculating the tangents
    dcur = cell(Nn-1,3);
    
    for j = 1:(Nn-1)
        dcur{j,3} = (x{j+1} - x{j})/norm(x{j+1} - x{j});
    end
    
    
    
    %% Calculating the twist quaternion
    ptwist = cell(Nn-1,1);
    
    for j = 1:(Nn-1)
        r = cos(gamma(j)/2);
        v = sin(gamma(j)/2)*dstar{j,3};
        ptwist{j} = qgen(r,v(1),v(2),v(3));
    end
    
    %% Calculating parallel transport quaternion
    ppar = cell(Nn-1,1);
    
    for j = 1:(Nn-1)
        
        r = sqrt((1+dstar{j,3}'*dcur{j,3})/2);
        v = (1/sqrt(2))*(1/sqrt(1 + dstar{j,3}'*dcur{j,3})) * cross(dstar{j,3},dcur{j,3});
        ppar{j} = qgen(r,v(1),v(2),v(3));
    end
    %% Calculating the product of quaternions
    pi = cell(Nn-1,1);
    
    for j = 1:(Nn-1)
        pi{j} = qmult(ppar{j},ptwist{j});
    end
    %% Calculating the directors
    
    for j = 1:(Nn-1)
        dcur{j,1} = quat2vec(qmult(pi{j},vec2quat(dstar{j,1}),qconj(pi{j})));
        dcur{j,2} = quat2vec(qmult(pi{j},vec2quat(dstar{j,2}),qconj(pi{j})));
        dcur{j,3} = (x{j+1} - x{j})/norm(x{j+1} - x{j});
    end
    
    
   [omega,eta] = CalcOmegaEta(Xextender(Xsolve,dS),Sys);
   
    edgevecsp = cell(Nn,1); 
    edgevecsm = cell(Nn,1);
    for j = 2:(Nn-1)
        edgevecsp{j} = w*(dcur{j,1} + eta(j)*dcur{j,3});
    end
    edgevecsp{1} = edgevecsp{2}; edgevecsp{end} = edgevecsp{end-1};
    for j = 1:Nn
        edgevecsm{j} = -edgevecsp{j};        
    end
    
    
    
    %%
    %tocstep = toc;
    if ploton
        plot3(Xvals(:,1),Xvals(:,2),Xvals(:,3));
        hold all
        plot3(Xvals(:,1),Xvals(:,2),Xvals(:,3),'b.');
        plot3(Xseg(:,1),Xseg(:,2),Xseg(:,3),'r.');
        
        for j = 1:(Nn-1)
            plot3([Xseg(j,1);Xseg(j,1) + (1/2)*sc*dcur{j,1}(1)],[Xseg(j,2);Xseg(j,2) + (1/2)*sc*dcur{j,1}(2)],[Xseg(j,3);Xseg(j,3) + (1/2)*sc*dcur{j,1}(3)],'r-')
            plot3([Xseg(j,1);Xseg(j,1) + (1/2)*sc*dcur{j,2}(1)],[Xseg(j,2);Xseg(j,2) + (1/2)*sc*dcur{j,2}(2)],[Xseg(j,3);Xseg(j,3) + (1/2)*sc*dcur{j,2}(3)],'g-')
        end
        
        
        for j = 1:Nn
            if j == 1
                plotCircle3D2(x{j},dcur{j,1},dcur{j,2},(1/2)*sc)
            elseif j == Nn
                plotCircle3D2(x{j},dcur{j-1,1},dcur{j-1,2},(1/2)*sc)
            else
                plotCircle3D2(x{j},(dcur{j-1,1} + dcur{j,1})/2,(dcur{j-1,2} + dcur{j,2})/2,(1/2)*sc)
            end
        end
        
        
        hold off
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %%
        title([num2str(i) ', Degrees of Freedom: ' num2str(max(size(Xsolve))) ', TocStep: ' num2str(tocstep) ', E = ' num2str(ff.objective(Xsolve))])
        axis equal
        drawnow
    end
    if savedata
        xv = zeros(Nn,1);
        yv = zeros(Nn,1);
        zv = zeros(Nn,1);
        d11 = zeros(Nn,1);
        d12 = zeros(Nn,1);
        d13 = zeros(Nn,1);
        for ii = 1:Nn
            
            xv(ii) = x{ii}(1);
            yv(ii) = x{ii}(2);
            zv(ii) = x{ii}(3);
            if ii < Nn
                d11(ii) = dcur{ii,1}(1);
                d12(ii) = dcur{ii,1}(2);
                d13(ii) = dcur{ii,1}(3);
            end
        end
        MM = [xv,yv,zv,d11,d12,d13];
        TT = array2table(MM,'VariableNames',{'xcoord','ycoord','zcoord','d11','d12','d13'});
        
        %writetable(TT,['data/data.' num2str(i) '.csv']);
        %vtkwrite(['data/data.' num2str(i) '.vtk'], 'polydata', 'lines', xv, yv, zv)
        datapad = @(X) [X(2);X(2:(end-1));X(end-1)];
        writevtkline(['data/data.' num2str(i) '.vtk'],xv,yv,zv,'omega', datapad(omega)/dS ,'edgevecsp',edgevecsp)
    end
    dstar = dcur;
    %tocstep
end

%% Axial strain
axstr = zeros(Nn-1,1);
for i = 1:(Nn-1)
    axstr(i) = norm(x{i+1} - x{i});
end
figure()
plot(Si(1:(end-1)),axstr-dS)
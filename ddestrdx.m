function ddestrout = ddestrdx(ti,epsvali,dS)


ddestrout = zeros(6,6);
% for i = 1:3
%     for j = 1:3
%         %dxidxi
%         d1xi = zeros(3,1); d1xi(i) = 1;
%         d2xi = zeros(3,1); d2xi(j) = 1;
%         d2ti = -(eye(3) - ti*ti')*d2xi;
%         ddestrout(i,j) = (ti'*d2xi)*(ti'*d1xi) - (epsvali - dS)*d2ti'*(d1xi);
%         %dxidxip
%         d1xi = zeros(3,1); d1xi(i) = 1;
%         d2xip = zeros(3,1); d2xip(j) = 1;
%         d2ti = (eye(3) - ti*ti')*d2xip;
%         ddestrout(i,j+3) = -(ti'*d2xip)*(ti'*d1xi) - (epsvali - dS)*d2ti'*(d1xi);
%         ddestrout(j+3,i) = ddestrout(i,j+3);
%         
%         
%         %dxipdxip
%         d1xip = zeros(3,1); d1xip(i) = 1;
%         d2xip = zeros(3,1); d2xip(j) = 1;
%         d2ti = (eye(3) - ti*ti')*d2xip;
%         ddestrout(i+3,j+3) = (ti'*d2xip)*(ti'*d1xip) + (epsvali - dS)*d2ti'*d1xip;
%         
%         
%     end
% end

%dxidxi
ddestrout(1:3,1:3) = ti*ti' - (epsvali - dS)*(-(eye(3) - ti*ti'))/epsvali;
ddestrout(1:3,4:6) = -ti*ti' - (epsvali - dS)*(eye(3) - ti*ti')/epsvali;
ddestrout(4:6,1:3) = ddestrout(1:3,4:6);
ddestrout(4:6,4:6) = ti*ti' + (epsvali - dS)*(eye(3) - ti*ti')/epsvali;
end
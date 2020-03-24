function [omega,eta] = kappa2omegaeta(kappa1,kappa3,dS)
% N = max(size(kappa1));
% omega = kappa1;
% eta = zeros(size(omega));
% 
% etadefined = (omega/dS)> 0.000001;
% for i = 1:N
%    if etadefined(i)
%        eta(i) = kappa3(i)/kappa1(i);
%    end
% end
% etadefined(1) = 1; etadefined(end) = 1;
% 
% i = 2;
% while i < N
%    if etadefined(i) == 0
%        x2 = findnext(i);
%        x1 = i-1;
%        
%        y2 = eta(x2);
%        y1 = eta(x1);
%        
%        M = (y2-y1)/(x2-x1);
%        interphan = @(k) y1 + M*(k-x1);
%        
%        for j = i:(x2-1)
%            eta(j) = interphan(j);
%        end
%        i = x2;
%        
%    end
%    i = i+1;
% end
% 
%     function iout = findnext(j)
%         iout = 0;
%         while iout == 0
%             if etadefined(j)
%                 iout = j;
%             else 
%                 j=j+1;
%             end
%             
%         end
%     end
N = max(size(kappa1));
omega = kappa1;
eta = zeros(size(kappa1));
for i = 1:N
    if abs(omega(i))>0
        eta(i) = kappa3(i)/kappa1(i);
    end
end

end


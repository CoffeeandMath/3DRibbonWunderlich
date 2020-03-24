function ddEbendiout = ddEbendi(~,Sys)
Nn = Sys.Nn;
D = Sys.D;


ddEbendiout = zeros(3,3,Nn);
for i = 2:(Nn-1)
    ddEbendiout(:,:,i) = diag(D(i,:));
end

end


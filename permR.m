function permRout = permR(R)
permRout = zeros(3,1);
permRout(1) = R(2,3) - R(3,2);
permRout(2) = R(3,1) - R(1,3);
permRout(3) = R(1,2) - R(2,1);
end


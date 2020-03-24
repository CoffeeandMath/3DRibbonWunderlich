function qlogout = qlog(q)

[a,b,c,d] = parts(q);

qlog1 = quatlog(quatnormalize([a b c d]));
qlogout = quaternion(qlog1(1),qlog1(2),qlog1(3),qlog1(4));

end


function err = testGeodesics(k ,a ,b ,c ,d ,alpha ,beta ,z)

z1 = z;
[Phi ,G] = Geodesics(a ,b ,c ,d ,alpha ,beta ,z);

for i = 1:5
    z1(k) = z(k)+10^(-i);
    [phi ,g] = Geodesics(a ,b ,c ,d ,alpha ,beta ,z1);
    err(i,1) = G(k) - (phi-Phi)/10^(-i);
end




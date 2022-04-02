
% dataA      = readmatrix("Aguilar.csv");
wavelength = (0.3:0.001:1);
% ZnOindex   = spline(dataA(1:581,1),dataA(1:581,2),wavelength);
IPAindex   = IPA_f(wavelength);


function f = SiO2_f(lambda) %Malitson 1965: Fused silica; n 0.21-6.7 µm
f1 =  (0.6961663.*(lambda.^2))./((lambda.^2)-(0.0684043^2));
f2 =  (0.4079426.*(lambda.^2))./((lambda.^2)-(0.1162414^2));
f3 =  (0.8974794.*(lambda.^2))./((lambda.^2)-(9.896161^2));
f  = sqrt(f1+f2+f3+1);
end

function f = IPA_f(lambda) %Sani and Dell'Oro 2016: iso-propanol;Sellmier; n,k 0.185-2.8 µm
f1 =  (0.0107.*(lambda.^2))./((lambda.^2)-(8.88^2));
f2 =  (0.8702.*(lambda.^2))./((lambda.^2)-(0.01036^2));
f  = sqrt(f1+f2+1);
end

function f = ZnO_f(lambda) %ZnOthinflim_2005; Cauchy; Unknown range
f1 =  0.059./(lambda.^2);
f  = f1+.781;
end


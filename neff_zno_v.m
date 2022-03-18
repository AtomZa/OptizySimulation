clear all
close all
clc
%-----------------------------------------------------------------------------%
nmed = 1.0;
lambdavis = 0.6:0.001:0.7;
lambdair = 1.5:0.001:1.6;
lambdamainvis = 0.633;
volfact = 0.0:0.1:1.0; %Change volume fraction of interest here

nvis = zeros(size(lambdavis));
for nn = 1:length(lambdavis)
    nvis(nn) = 1.781 + (0.059/(lambdavis(nn)^2));
end

n633 = 1.781 + (0.059/(lambdamainvis^2));

nir = zeros(size(lambdair));
for nn = 1:length(lambdair)
    nir(nn) = 1.781 + (0.059/(lambdair(nn)^2));
end

neffvis = zeros(length(volfact),length(lambdavis));
for nn = 1:length(volfact)
    for mm = 1:length(lambdavis)
        neffvis(nn,mm) = (volfact(nn)*nvis(mm)) + ((1 - volfact(nn))*nmed);
    end
end

neffir = zeros(length(volfact),length(lambdair));
for nn = 1:length(volfact)
    for mm = 1:length(lambdair)
        neffir(nn,mm) = (volfact(nn)*nir(mm)) + ((1 - volfact(nn))*nmed);
    end
end

neff633 = zeros(length(volfact),1);
for nn = 1:length(volfact)
    neff633(nn,1) = (volfact(nn)*n633) + ((1 - volfact(nn))*nmed);
end
save('neff_zno_data_10%.mat','lambdair','lambdavis','neff633','neffir','neffvis','volfact','nvis','n633','nir')
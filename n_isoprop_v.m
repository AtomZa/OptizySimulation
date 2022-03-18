clear all
close all
clc
%-----------------------------------------------------------------------------%
nmed = 1.00;
lambdavis = 0.6:0.001:0.7;
lambdair = 1.5:0.001:1.6;
lambdamainvis = 0.633;
volfact = 0:0.1:1; %Change volume fraction of interest here

nvis = zeros(size(lambdavis));
for nn = 1:length(lambdavis)
    nvis(nn) = sqrt(1 + (0.0107*(lambdavis(nn)^2)/((lambdavis(nn)^2) - 8.88)) + (0.8702*(lambdavis(nn)^2)/((lambdavis(nn)^2) - 0.01036)));
end

n633 = sqrt(1 + (0.0107*(lambdamainvis^2)/((lambdamainvis^2) - 8.88)) + (0.8702*(lambdamainvis^2)/((lambdamainvis^2) - 0.01036)));

nir = zeros(size(lambdair));
for nn = 1:length(lambdair)
    nir(nn) = sqrt(1 + (0.0107*(lambdair(nn)^2)/((lambdair(nn)^2) - 8.88)) + (0.8702*(lambdair(nn)^2)/((lambdair(nn)^2) - 0.01036)));
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
save('neff_isoprop_data_10%.mat','lambdair','lambdavis','neff633','neffir','neffvis','volfact','nvis','n633','nir')
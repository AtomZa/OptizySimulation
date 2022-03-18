clear all
close all
clc

nmed = 1.0;
lambdavis = 0.4:0.001:0.8;
lambdair = 1.54:0.001:1.56;
% lambdamainvis = 0.633;
volfact = 0.003:0.1:0.003; %Change volume fraction of interest here


n633 = 0.18344 + (3.4332*1i); %data from refractiveindex.info

nir = zeros(size(lambdair));
for nn = 1:length(lambdair)
    tempn = 0.0445*(exp(1.5954*lambdair(nn)));
    tempk = (7.7741*lambdair(nn)) - 1.2975;
    nir(nn) = tempn + (tempk*1i);
end

nvis = zeros(size(lambdavis));
for nn = 1:length(lambdair)
    tempn = 0.0445*(exp(1.5954*lambdair(nn)));
    tempk = (7.7741*lambdair(nn)) - 1.2975;
    nir(nn) = tempn + (tempk*1i);
end

% neffvis = zeros(length(volfact),length(lambdavis));
% for nn = 1:length(volfact)
%     for mm = 1:length(lambdavis)
%         neffvis(nn,mm) = (volfact(nn)*nvis(mm)) + ((1 - volfact(nn))*nmed);
%     end
% end

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
save('neff_gold_data.mat','lambdair','neffir','volfact','nir','n633')
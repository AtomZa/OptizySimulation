clear all
close all
clc
%-----------------------------------------------------------------------------%
global ncoresmf ncladdingsmf asmf ncoremmf ncladdingmmf ammf zmmf epscoresmf epscladdingsmf epscoremmf epscladdingmmf k0 neffarrsmf lambda nmedium aair
ncoresmf = 1.4504; %core refractive index of smf
ncladdingsmf = 1.447; %cladding refractive index of smf
asmf = 4.1; %radius of smf
ncoremmf = 1.443; %refractive index of no core mmf
nmedium = 1.04641231840654; %medium refractive index
ncladdingmmf = nmedium; %refractive index of medium surrounding no core mmf
ammf = 62.5; %radius of no core mmf
zmmf = 57713; %length of mmf segment
epscoresmf = ncoresmf^2; %core pzermittivity of smf
epscladdingsmf = ncladdingsmf^2; %cladding permittivity of smf
epscoremmf = ncoremmf^2; %permittivity of no core mmf
epscladdingmmf = ncladdingmmf^2; %permittivity of medium surrounding no core mmf
% lambda = 1.55; %wavelength
lambda = 0.633;
k0 = 2*pi/lambda; %wave number
aair = 87.5; %radius of medium
%----------------------finding mode in smf--------------------------------%
increment = 0.0001;
ifinal = ((ncoresmf - ncladdingsmf)/increment);
neffarrsmf = zeros(1,round(ifinal)); %empty array of effective index of smf to be assigned in the next loop
for i = 0:ifinal
     try 
        chk = fzero(@characeqnsmf,[(ncladdingsmf + (i*increment)),(ncladdingsmf + ((i+1)*increment))]);
     catch
        neffarrsmf(i+1) = 0;
        continue
     end
        neffarrsmf(i+1) = fzero(@characeqnsmf,[(ncladdingsmf + (i*increment)),(ncladdingsmf + ((i+1)*increment))]);
end

neffarrsmf = neffarrsmf(neffarrsmf ~= 0);
neffarrsmf = neffarrsmf(length(neffarrsmf));

rarrsmf = -asmf:0.01:asmf; %array of r from -asmf to asmf
usmf = findusmf(neffarrsmf); %u of smf
wsmf = findwsmf(neffarrsmf); %w of smf
rarrcladding1 = -ammf:0.01:(-asmf-0.01); %array of r from -ammf to -asmf
rarrcladding2 = (asmf + 0.01):0.01:ammf; %array of r from asmf to ammf
rarrair1 = -aair:0.01:(-ammf-0.01); %array of r from -aair to -ammf
rarrair2 = (ammf + 0.01):0.01:aair; %array of r frin ammf to aair
rarrall = [rarrair1 rarrcladding1 rarrsmf rarrcladding2 rarrair2]; %array of r from -aair to aair

ecoresmf = besselj(0,(usmf.*rarrsmf./asmf)); %field where -asmf < r < asmf
c = (besselj(0,usmf))/(besselk(0,wsmf)); %scaling coefficient
ecladdingsmf2 = c.*(besselk(0,(wsmf.*rarrcladding2./asmf))); %field where asmf < r < ammf
ecladdingsmf1 = zeros(size(ecladdingsmf2)); %field where -ammf < r < -asmf
for nn = 1:length(ecladdingsmf1)
    ecladdingsmf1(nn) = ecladdingsmf2(length(ecladdingsmf2) - (nn - 1));
end
eair1 = zeros(size(rarrair1)); %field where -aair < r < -ammf
eair2 = zeros(size(rarrair2)); %field where ammf < r < aair
eair1(:) = ecladdingsmf1(1);
eair2(:) = ecladdingsmf2(length(ecladdingsmf2));
eall = [eair1 ecladdingsmf1 ecoresmf ecladdingsmf2 eair2]; %field where -aair < r < aair
iall = zeros(size(eall)); %field intensity where -aair < r < aair
for nn = 1:length(eall)
    iall(nn) = (abs(eall(nn)))^2;
end
figure()
plot(rarrall,eall)
%----------------------finding mode in mmf--------------------------------%
ifinal = ((ncoremmf - ncladdingmmf)/increment);
neffarrmmf = zeros(1,round(ifinal)); %empty array of effective indices of mmf to be assigned in the next loop
for ii = 0:ifinal
    try 
        chk = fzero(@characeqnmmf,[(ncladdingmmf + (ii*increment)),(ncladdingmmf + ((ii+1)*increment))]);
    catch
         neffarrmmf(ii+1) = 0;
         continue
    end
         neffarrmmf(ii+1) = fzero(@characeqnmmf,[(ncladdingmmf + (ii*increment)),(ncladdingmmf + ((ii+1)*increment))]);
end
neffarrmmf = neffarrmmf(neffarrmmf ~= 0);
for p = 1:(length(neffarrmmf))/2
    temp = neffarrmmf(p);
    neffarrmmf(p) = neffarrmmf(length(neffarrmmf) - (p - 1));
    neffarrmmf(length(neffarrmmf) - (p - 1)) = temp;
end
figure()
plot(neffarrmmf,'.')
%------------------------finding u and w for each mode in mmf----------------------------------%
uarrmmf = zeros(size(neffarrmmf)); %u of each mode in mmf
warrmmf = zeros(size(neffarrmmf)); %w of each mode in mmf
carr = zeros(size(neffarrmmf)); %scaling coefficient
for j = 1:length(neffarrmmf)
    uarrmmf(j) = findummf(neffarrmmf(j));
    warrmmf(j) = findwmmf(neffarrmmf(j));
    carr(j) = (besselj(0,uarrmmf(j)))/(besselk(0,warrmmf(j),1));
end
%------------------------finding field and its intensity for each mode in mmf----------------------------------%
earrallmmf = zeros(length(neffarrmmf),length(rarrall)); %field of each mode in mmf

for k = 1:length(neffarrmmf)
    for m = 1:length(rarrall)
        if (rarrall(m) > ammf)
            earrallmmf(k,m) = carr(k)*(besselk(0,(warrmmf(k)*rarrall(m)/ammf),1));
        elseif (rarrall(m) < -ammf)
            earrallmmf(k,m) = carr(k)*(besselk(0,(warrmmf(k)*rarrall(length(rarrall) - (m - 1))/ammf),1));
        else
            earrallmmf(k,m) = besselj(0,(uarrmmf(k)*rarrall(m)/ammf));
        end
    end
end

iarrallmmf = zeros(size(earrallmmf)); %intensity of each mode in mmf
for nn = 1:length(neffarrmmf)
    for mm = 1:length(rarrall)
        iarrallmmf(nn,mm) = (abs(earrallmmf(nn,mm)))^2;
    end
end
%-----------------------finding variable pj and aj for each mode in mmf-----------------------------------------%
pj = zeros(size(neffarrmmf)); %empty array of parameter pj 
aj = zeros(size(neffarrmmf)); %empty array of parameter aj
for nn = 1:length(neffarrmmf)
    funpj1 = @(r) ((abs(besselj(0,(uarrmmf(nn).*r/ammf)))).^2).*r; %function of pj where r is between 0 and ammf
    funpj2 = @(r) ((abs(carr(nn).*besselk(0,(warrmmf(nn).*r/ammf)))).^2).*r; %function of pj where r is greater than ammf
    pj(nn) = ((integral(funpj1,0,ammf)) + (integral(funpj2,ammf,Inf))).*2.*pi;
end

integ1 = zeros(size(neffarrmmf)); %empty arrays for storing integrals in the next loop
integ2 = zeros(size(neffarrmmf));
integ3 = zeros(size(neffarrmmf));
for nn = 1:length(neffarrmmf)
%     funaj1 = @(r) (cross((besselj(0,(usmf.*r./asmf))),conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;
%     funaj2 = @(r) (cross((besselk(0,(usmf.*r./asmf),1)),conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;
%     funaj3 = @(r) (cross(eair2(length(eair2)),conj(besselk(0,(uarrmmf(nn).*r/ammf),1)))).*r;
% error with cross product so scalar multiplication is applied instead %
    funaj1 = @(r) ((besselj(0,(usmf.*r./asmf))).*(conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r; %function where r is between 0 and asmf
    funaj2 = @(r) ((c.*besselk(0,(wsmf.*r./asmf))).*(conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r; %function where r is between asmf and ammf
    funaj3 = @(r) ((eair2(length(eair2))).*(conj(carr(nn).*besselk(0,(warrmmf(nn).*r/ammf))))).*r; %function where r is greater than ammf
    integ1(nn) = (integral(funaj1,0,asmf)); %integrate from 0 to asmf
    integ2(nn) = (integral(funaj2,asmf,ammf)); %integrate from asmf to ammf
    integ3(nn) = (integral(funaj3,ammf,Inf)); %integrate from ammf to inf
    aj(nn) = (integ1(nn) + integ2(nn) + integ3(nn)).*2.*pi./pj(nn);
end
%-----------------------finding efield of mmf with interference---------------------------------------%
betammf = neffarrmmf.*k0; %array of propagation constants of each mode
% zarr = 57000:10:60000; %array of propagation length
zarr = 50000:100:70000;
%change zarr if the code take too long to run
einter = zeros(length(zarr),length(rarrall)); %array of field of multimode interference
einter0 = zeros(1,length(rarrall)); %array of field of multimode interference at z = 0
for nn = 1:length(zarr)
    for mm = 1:length(rarrall)
        for pp = 1:length(neffarrmmf)
            einter(nn,mm) = einter(nn,mm) + (aj(pp)*earrallmmf(pp,mm)*exp(1j*betammf(pp)*zarr(nn)));
        end
    end
end
for mm = 1:length(rarrall)
    for pp = 1:length(neffarrmmf)
        einter0(mm) = einter0(mm) + (aj(pp)*earrallmmf(pp,mm)*exp(1j*betammf(pp)*0));
    end
end
iinter = (abs(einter)).^2; %array of field intensity of multimode interference
iinter0 = (abs(einter0)).^2; %array of field intensity of multimode interference at z = 0
figure()
contourf(rarrall,zarr,iinter,'LineStyle','none','LevelList',0:0.01:1)
hold on
colormap(jet)

%-----------------------finding coupling efficiency of mmf with interference----------------------------------%
funescore = @(r) ((abs(besselj(0,(usmf.*r./asmf)))).^2).*r;
funesclad = @(r) ((abs(c*besselk(0,(wsmf.*r./asmf)))).^2).*r;
funesair = @(r) ((ecladdingsmf2(length(ecladdingsmf2))).^2).*r;
% ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + integral(funesair,ammf,Inf)).*(2*pi);
ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + 0).*(2*pi); %assume third integral is zero
ajbar = aj.*(sqrt(pj./ps)); %aj bar parameters of each mode

etaarr = zeros(size(zarr));
for nn = 1:length(zarr)
    for mm = 1:length(neffarrmmf)
        beth = betammf(mm);
        ahb = ajbar(mm);
        for pp = 1:length(neffarrmmf)
            ajb = ajbar(pp);
            betj = betammf(pp);
            etaarr(nn) = etaarr(nn) + ((ajb^2)*((conj(ahb))^2)*exp(1j*(betj - beth)*zarr(nn)));
        end
    end
end
figure()
plot(zarr,10.*log(etaarr)./log(10))
eta0 = 0;
for mm = 1:length(neffarrmmf)
    beth = betammf(mm);
    ahb = ajbar(mm);
    for pp = 1:length(neffarrmmf)
        ajb = ajbar(pp);
        betj = betammf(pp);
        eta0 = eta0 + ((ajb^2)*((conj(ahb))^2)*exp(1j*(betj - beth)*0));
    end
end
% figure()
% plot(zarr,10.*log(etaarr)./log(10))
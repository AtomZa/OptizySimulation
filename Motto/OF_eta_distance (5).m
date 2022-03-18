clearvars
close all
clc
%-----------------------------------------------------------------------------%
global ncoresmf ncladdingsmf asmf ncoremmf ncladdingmmf ammf epscoresmf epscladdingsmf epscoremmf epscladdingmmf k0 neffarrsmf lambda nmedium aair
ncoresmf = 1.4504;                     %core refractive index of smf
ncladdingsmf = 1.447;                 %cladding refractive index of smf
asmf = 4.1;                                 %radius of smf
ncoremmf = 1.443;                      %refractive index of no core mmf
nmedium = 1.04641231840654;   %medium refractive index
ncladdingmmf = nmedium;            %refractive index of medium surrounding no core mmf
ammf = 62.5;                                %radius of no core mmf
%zmmf = 100000;                              %length of mmf segment
epscoresmf = ncoresmf^2;             %core pzermittivity of smf
epscladdingsmf = ncladdingsmf^2; %cladding permittivity of smf
epscoremmf = ncoremmf^2;           %permittivity of no core mmf
epscladdingmmf = ncladdingmmf^2;%permittivity of medium surrounding no core mmf
lambda = 0.633;
k0 = 2*pi/lambda;                         %wave number
aair = 87.5;                                   %radius of medium
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
% figure()
% plot(rarrall,eall)
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
% figure()
% plot(neffarrmmf,'.')
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

% error with cross product so scalar multiplication is applied instead %
    funaj1 = @(r) ((besselj(0,(usmf.*r./asmf))).*(conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;       %function where r is between 0 and asmf
    funaj2 = @(r) ((c.*besselk(0,(wsmf.*r./asmf))).*(conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;  %function where r is between asmf and ammf
    funaj3 = @(r) ((eair2(length(eair2))).*(conj(carr(nn).*besselk(0,(warrmmf(nn).*r/ammf))))).*r; %function where r is greater than ammf
    integ1(nn) = (integral(funaj1,0,asmf));       %integrate from 0 to asmf
    integ2(nn) = (integral(funaj2,asmf,ammf)); %integrate from asmf to ammf
    integ3(nn) = (integral(funaj3,ammf,Inf));     %integrate from ammf to inf
    aj(nn) = (integ1(nn) + integ2(nn) + integ3(nn)).*2.*pi./pj(nn);
end
%-----------------------finding efield of mmf with interference---------------------------------------%

% ---required element    --- %
%            betammf
%            earrallmmf
%            zarr(arbitary)


betammf = gpuArray(neffarrmmf.*k0); %array of propagation constants of each mode
zarr = gpuArray(50000:50:145000);
rres = (7000:8750);
aj = gpuArray(aj);
earrallmmf = gpuArray(earrallmmf(:,rres));

disp 'mmf interferer'
tic

nn   = int16(1:length(zarr));
mm = int16(1:length(rarrall(rres)));
   
p1 = (earrallmmf(:,mm)).*aj.';
p2 = exp(1j*(betammf.').*(zarr(nn))); 

einter = (p2.')*p1;
iinter = gather(abs(einter)).^2; %array of field intensity of multimode interference
einter = gather(einter.');
toc


% disp 'countour time'
% tic
%  figure()
% contourf(rarrall(5751:11750),zarr,iinter,'LineStyle','none','LevelList',0:0.005:1);
% hold on
% colormap(jet)
% toc

 figure()
 rescaledMatrix = uint8(rescale(iinter, 0, 255));
 doubleRe = flip([rescaledMatrix flip(rescaledMatrix,2)]);
 imshow(doubleRe,jet);


%-----------------------finding coupling efficiency of mmf with interference----------------------------------%

% ---required element    --- %
%            aj
%            pj
%            zarr(arbitary)
%--- required variable ---- %
%            usmf
%            wsmf

funescore = @(r) ((abs(besselj(0,(usmf.*r./asmf)))).^2).*r;
funesclad = @(r) ((abs(c*besselk(0,(wsmf.*r./asmf)))).^2).*r;
ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf)).*(2*pi); %assume third integral is zero
ajbar = aj.*(sqrt(pj./ps));
ajbar = ajbar(1:25);
betammf1 = betammf(1:25);
disp('my new method (coupling eff)')
tic

nn = 1:length(zarr);

P100= reshape(((betammf1 - betammf1.')*1j).',[],1);
P110 = P100.*zarr(nn);
P220= reshape(log((ajbar.^2).*((conj(ajbar)).^2).'),1,[]);
P330 = sum(exp(P220+P110.'),2);

c1oup =  abs(reshape(sum(P330,2),[],1)).';
toc

figure()
plot(zarr,10.*log(c1oup)./log(10))

%first one%
% ajbar = gather(aj.*(sqrt(pj./ps))); %aj bar parameters of each mode
% ajbar = ajbar(1:20);
% betammf = gather(betammf(1:20));
% zarr = gather (zarr);

% disp('my first method')
% tic
% 
% metaarr = zeros(length(zarr));
% for nn = 1:length(zarr)
%     P11= exp(1j*(betammf - betammf.').*(zarr(nn))).';
%     P22= ((ajbar.^2).*((conj(ajbar)).^2).');
%     metaarr(nn) = sum(P22.*P11,"all");
% end
% 
% coup =  abs(reshape(sum(metaarr,2),[],1)).';
% toc
% 
% figure()
% plot(zarr,10.*log(coup)./log(10))

% disp('reduction dimension')
% tic
% etaarr = zeros(length(zarr));
% for nn = 1:length(zarr)
%     for mm = 1:25  % only 25~50 element are matter (err = -5*10^-8 {1:25})
%         for pp = 1:25                              
%             etaarr(nn) = etaarr(nn) + (((ajbar(pp)^2)*(conj(ajbar(mm)))^2)*exp(1j*(betammf(pp) - betammf(mm))*zarr(nn)));
%         end 
%     end
% end
% e2taarr = real(reshape(sum(etaarr,2),[],1).');
% toc
% 
% figure()
% plot(zarr,10.*log(e2taarr)./log(10))
% 
% disp('old one')
% tic
% etaarr = zeros(length(zarr));
% for nn = 1:length(zarr)
%     for mm = 1:length(neffarrmmf)
%         for pp = 1:length(neffarrmmf)                              
%             etaarr(nn) = etaarr(nn) + (((ajbar(pp)^2)*(conj(ajbar(mm)))^2)*exp(1j*(betammf(pp) - betammf(mm))*zarr(nn)));
%         end 
%     end
% end
% etaarr = real(reshape(sum(etaarr,2),[],1).');
% toc
% 
% figure()
% plot(zarr,10.*log(etaarr)./log(10))
% 
% diff = etaarr-e2taarr;
% diff = sum(diff,"all")/length(diff);
% disp(diff)


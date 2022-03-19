clear all
close all
clc
% save('main_new_1500_1600_n12.mat','etaarr','zmmf','nmedium','aair','lambdaarr')
%-----------------------------------------------------------------------------%
global ncoresmf ncladdingsmf asmf ncoremmf ncladdingmmf ammf zmmf epscoresmf epscladdingsmf epscoremmf epscladdingmmf k0 neffarrsmf lambda nmedium aair
zno_data = load('neff_zno_data.mat');
isoprop_data = load('neff_isoprop_data.mat');
nzno = zno_data.nir;
nsolution = isoprop_data.neffir;
volfact = zno_data.volfact;
lambdaarr = 1.50:0.001:1.60;
nmedium = cell(length(volfact),length(volfact)); %medium refractive index
for mm = 1:length(volfact)
    for qq = 1:length(volfact)
        tempnmedium = zeros(1,length(lambdaarr));
        for nn = 1:length(lambdaarr)
            tempnmedium(1,nn) = (volfact(qq)*nzno(nn)) + ((1 - volfact(qq))*nsolution(mm,nn));
        end
    nmedium{mm,qq} = tempnmedium;    
    end
end
ncoresmf = 1.4504; %core refractive index of smf
ncladdingsmf = 1.447; %cladding refractive index of smf
asmf = 4.1; %radius of smf
ncoremmf = 1.443; %refractive index of no core mmf
ncladdingmmf = nmedium; %refractive index of medium surrounding no core mmf
ammf = 62.5; %radius of no core mmf
zmmf = 58540; %length of mmf segment
epscoresmf = ncoresmf^2; %core permittivity of smf
epscladdingsmf = ncladdingsmf^2; %cladding permittivity of smf
epscoremmf = ncoremmf^2; %permittivity of no core mmf
% epscladdingmmf = ncladdingmmf.^2; %permittivity of medium surrounding no core mmf
% k0 = 2*pi/lambda; %wave number
aair = 87.5; %radius of medium
% lambdaarr = 0.4:0.001:0.5;
k0arr = 2.*pi./lambdaarr;

%----------------------finding mode in smf--------------------------------%

increment = 0.0001;
ifinal = ((ncoresmf - ncladdingsmf)/increment);
neffarrsmf = zeros(length(lambdaarr),round(ifinal));
for nn = 1:length(lambdaarr)
    k0 = k0arr(nn);
    for i = 0:ifinal
        try 
            chk = fzero(@characeqnsmf,[(ncladdingsmf + (i*increment)),(ncladdingsmf + ((i+1)*increment))]);
        catch
            neffarrsmf(nn,i+1) = 0;
            continue
        end
        neffarrsmf(nn,i+1) = fzero(@characeqnsmf,[(ncladdingsmf + (i*increment)),(ncladdingsmf + ((i+1)*increment))]);
    end
end
temp = zeros(size(lambdaarr));
for nn = 1:length(lambdaarr)
    for mm = 1:length(neffarrsmf(1,:))
        if neffarrsmf(nn,mm) ~= 0
            temp(nn) = neffarrsmf(nn,mm);
        end
    end
end
neffarrsmf = temp;
% plot(lambdaarr,neffarrsmf)

usmf = zeros(size(neffarrsmf));
wsmf = zeros(size(neffarrsmf));
csmf = zeros(size(neffarrsmf));
for nn = 1:length(lambdaarr)
    k0 = k0arr(nn);
    usmf(nn) = findusmf(neffarrsmf(nn));
    wsmf(nn) = findwsmf(neffarrsmf(nn));
    csmf(nn) = (besselj(0,usmf(nn)))/(besselk(0,wsmf(nn)));
end

rarrsmf = -asmf:0.01:asmf; %array of r from -asmf to asmf
rarrcladding1 = -ammf:0.01:(-asmf-0.01); %array of r from -ammf to -asmf
rarrcladding2 = (asmf + 0.01):0.01:ammf; %array of r from asmf to ammf
rarrair1 = -aair:0.01:(-ammf-0.01); %array of r from -aair to -ammf
rarrair2 = (ammf + 0.01):0.01:aair; %array of r frin ammf to aair
rarrall = [rarrair1 rarrcladding1 rarrsmf rarrcladding2 rarrair2]; %array of r from -aair to aair

eall = zeros(length(lambdaarr),length(rarrall));
for mm = 1:length(lambdaarr)
    ecoresmf = besselj(0,(usmf(mm).*rarrsmf./asmf)); %field where -asmf < r < asmf
    ecladdingsmf2 = csmf(mm)*besselk(0,(wsmf(mm).*rarrcladding2./asmf)); %field where asmf < r < ammf
    ecladdingsmf1 = zeros(size(ecladdingsmf2)); %field where -ammf < r < -asmf
    for nn = 1:length(ecladdingsmf1)
        ecladdingsmf1(nn) = ecladdingsmf2(length(ecladdingsmf2) - (nn - 1));
    end
    eair1 = zeros(size(rarrair1)); %field where -aair < r < -ammf
    eair2 = zeros(size(rarrair2)); %field where ammf < r < aair
    eair1(:) = ecladdingsmf1(1);
    eair2(:) = ecladdingsmf2(length(ecladdingsmf2));
    eall(mm,:) = [eair1 ecladdingsmf1 ecoresmf ecladdingsmf2 eair2]; %field where -aair < r < aair
end

iall = zeros(size(eall)); %field intensity where -aair < r < aair
for nn = 1:length(lambdaarr)
    for mm = 1:length(rarrall)
        iall(nn,mm) = (abs(eall(nn,mm)))^2;
    end
end

%----------------------finding mode in mmf--------------------------------%
neffcellmmf = cell(length(volfact),length(volfact));
for qq = 1:length(volfact)
    for ss = 1:length(volfact)
        tempcell = cell(1,length(lambdaarr));
        tempncladdingmmf = cell2mat(ncladdingmmf(qq,ss));
        if tempncladdingmmf(1) >= ncoremmf
            break
        end
        ifinal = zeros(1,length(lambdaarr));
        for rr = 1:length(ifinal)
            ifinal(rr) = ((ncoremmf - tempncladdingmmf(rr))/increment);
        end
        for nn = 1:length(lambdaarr)
            neffarrmmf = zeros(1,round(max(ifinal)));
            k0 = 2*pi/lambdaarr(nn);
            tempn = tempncladdingmmf(nn);
            epscladdingmmf = tempn^2;
            for ii = 0:ifinal(nn)
                try 
                    chk = fzero(@characeqnmmf,[(tempn + ((ii-1)*increment)),(tempn + (ii*increment))]);
                catch
                    neffarrmmf(ii+1) = 0;
                    continue
                end
                neffarrmmf(ii+1) = fzero(@characeqnmmf,[(tempn + ((ii-1)*increment)),(tempn + (ii*increment))]);
            end
            neffarrmmf = neffarrmmf(neffarrmmf ~= 0);
            for jj = 1:(length(neffarrmmf))/2
                temp = neffarrmmf(jj);
                neffarrmmf(jj) = neffarrmmf((length(neffarrmmf) - (jj - 1)));
                neffarrmmf((length(neffarrmmf) - (jj - 1))) = temp;
            end
            tempcell{1,nn} = neffarrmmf;
        end
        neffcellmmf{qq,ss} = tempcell;
    end
end
% counter = 0;
% for nn = 1:length(volfact)
%     if isempty(neffcellmmf{nn,1})
%         break
%     else
%         counter = counter + 1;
%     end
% end
% neffcellmmf = neffcellmmf(1:counter,1:length(lambdaarr));
% szneffcellmmf = size(neffcellmmf);

%-----------------------finding u and w of each mode in mmf-----------------------------------------%
wcellmmf = cell(size(neffcellmmf));
ucellmmf = cell(size(neffcellmmf));
ccellmmf = cell(size(neffcellmmf));
for qq = 1:length(neffcellmmf)
    for nn = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,nn})
            break
        end
        tempncladdingmmf = cell2mat(ncladdingmmf(qq,nn));
        tempneffcellmmf = neffcellmmf{qq,nn};
        tempucell = cell(1,length(lambdaarr));
        tempwcell = cell(1,length(lambdaarr));
        tempccell = cell(1,length(lambdaarr));
        for mm = 1:length(lambdaarr)
            k0 = 2*pi/lambdaarr(mm);
            tempn = tempncladdingmmf(mm);
            epscladdingmmf = tempn^2;
            temparr = cell2mat(tempneffcellmmf(1,mm));
            tempu = findummf(temparr);
            tempw = findwmmf(temparr);
            tempucell{1,mm} = tempu;
            tempwcell{1,mm} = tempw;
            tempccell{1,mm} = besselj(0,tempu)./besselk(0,tempw);
        end
        ucellmmf{qq,nn} = tempucell;
        wcellmmf{qq,nn} = tempwcell;
        ccellmmf{qq,nn} = tempccell;
    end
end

%-----------------------finding variable pj and aj for each mode in mmf-----------------------------------------%

pjcell = cell(size(neffcellmmf));
ajcell = cell(size(neffcellmmf));
for rr = 1:length(neffcellmmf)
    for qq = 1:length(neffcellmmf)
        if isempty(neffcellmmf{rr,qq})
            break
        end
        tempneffcellmmf = neffcellmmf{rr,qq};
        tempucellmmf = ucellmmf{rr,qq};
        tempwcellmmf = wcellmmf{rr,qq};
        tempccellmmf = ccellmmf{rr,qq};
        temppjcell = cell(1,length(lambdaarr));
        for mm = 1:length(lambdaarr)
            tempneff = cell2mat(tempneffcellmmf(1,mm));
            tempu = cell2mat(tempucellmmf(1,mm));
            tempw = cell2mat(tempwcellmmf(1,mm));
            tempc = cell2mat(tempccellmmf(1,mm));
            temppj = zeros(size(tempneff)); 
            for nn = 1:length(tempneff)
                funpj1 = @(r) ((abs(besselj(0,(tempu(nn).*r/ammf)))).^2).*r; %function of pj where r is between 0 and ammf
                funpj2 = @(r) ((abs(tempc(nn).*besselk(0,(tempw(nn).*r/ammf)))).^2).*r; %function of pj where r is greater than ammf
                temppj(nn) = ((integral(funpj1,0,ammf)) + (integral(funpj2,ammf,Inf))).*2.*pi;
            end
            temppjcell{1,mm} = temppj;
        end
        pjcell{rr,qq} = temppjcell;
    end
end
for rr = 1:length(neffcellmmf)
    for qq = 1:length(neffcellmmf)
        if isempty(neffcellmmf{rr,qq})
            break
        end
        tempneffcellmmf = neffcellmmf{rr,qq};
        tempucellmmf = ucellmmf{rr,qq};
        tempwcellmmf = wcellmmf{rr,qq};
        tempccellmmf = ccellmmf{rr,qq};
        temppjcell = pjcell{rr,qq};
        tempajcell = cell(1,length(lambdaarr));
        for mm = 1:length(lambdaarr)
            tempneff = cell2mat(tempneffcellmmf(1,mm));
            tempu = cell2mat(tempucellmmf(1,mm));
            tempw = cell2mat(tempwcellmmf(1,mm));
            tempc = cell2mat(tempccellmmf(1,mm));
            temppj = cell2mat(temppjcell(1,mm)); 
            integ1 = zeros(size(tempneff)); %empty arrays for storing integrals in the next loop
            integ2 = zeros(size(tempneff));
            integ3 = zeros(size(tempneff));
            tempaj = zeros(size(tempneff));
            for nn = 1:length(tempneff)
%     funaj1 = @(r) (cross((besselj(0,(usmf(mm).*r./asmf))),conj(besselj(0,(tempu(nn).*r/ammf))))).*r;
%     funaj2 = @(r) (cross((besselk(0,(usmf(mm).*r./asmf))),conj(besselj(0,(tempu(nn).*r/ammf))))).*r;
%     funaj3 = @(r) (cross(eall(mm,length(eall)),conj(besselk(0,(tempw(nn).*r/ammf))))).*r;
% error with cross product so scalar multiplication is applied instead %
                funaj1 = @(r) ((besselj(0,(usmf(mm).*r./asmf))).*(conj(besselj(0,(tempu(nn).*r/ammf))))).*r; %function where r is between 0 and asmf
                funaj2 = @(r) ((csmf(mm).*besselk(0,(wsmf(mm).*r./asmf))).*(conj(besselj(0,(tempu(nn).*r/ammf))))).*r; %function where r is between asmf and ammf
                funaj3 = @(r) (eall(mm,length(eall)).*(conj(tempc(nn).*besselk(0,(tempw(nn).*r/ammf))))).*r; %function where r is greater than ammf
                integ1(nn) = (integral(funaj1,0,asmf)); %integrate from 0 to asmf
                integ2(nn) = (integral(funaj2,asmf,ammf)); %integrate from asmf to ammf
                integ3(nn) = (integral(funaj3,ammf,Inf)); %integrate from ammf to inf
                tempaj(nn) = (integ1(nn) + integ2(nn) + integ3(nn)).*2.*pi./temppj(nn);
            end
            tempajcell{1,mm} = tempaj;
        end
        ajcell{rr,qq} = tempajcell;
    end
end
%-------------------------------------%
betacell = cell(size(neffcellmmf));
for qq = 1:length(neffcellmmf)
    for mm = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,mm})
            break
        end
        tempneffcellmmf = neffcellmmf{qq,mm};
        tempbetacell = cell(1,length(lambdaarr));
        for nn = 1:length(lambdaarr)
            k0 = 2*pi/lambdaarr(nn);
            tempneff = cell2mat(tempneffcellmmf(1,nn));
            tempbeta = tempneff.*k0;
            tempbetacell{1,nn} = tempbeta;
        end
        betacell{qq,mm} = tempbetacell;
    end
end

ps = zeros(size(neffarrsmf));
for nn = 1:length(lambdaarr)
    funescore = @(r) ((abs(besselj(0,(usmf(nn).*r./asmf)))).^2).*r;
    funesclad = @(r) ((abs(csmf(nn).*besselk(0,(wsmf(nn).*r./asmf)))).^2).*r;
%     funesair = @(r) ((eall(nn,length(eall))).^2).*r;
%     ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + integral(funesair,ammf,Inf)).*(2*pi);
    ps(nn) = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + 0).*(2*pi); %assume funesair is zero
end

ajbarcell = cell(size(neffcellmmf));
for rr = 1:length(neffcellmmf)
    for qq = 1:length(neffcellmmf)
        if isempty(neffcellmmf{rr,qq})
            break
        end
        temppjcell = pjcell{rr,qq};
        tempajcell = ajcell{rr,qq};
        tempajbarcell = cell(1,length(lambdaarr));
        for mm = 1:length(lambdaarr)
            tempps = ps(mm);
            temppj = cell2mat(temppjcell(1,mm));
            tempaj = cell2mat(tempajcell(1,mm));
            tempajbar = zeros(size(tempaj));
            for nn = 1:length(temppj)
                tempajbar(nn) = tempaj(nn)*(sqrt(temppj(nn)/tempps));
            end
            tempajbarcell{1,mm} = tempajbar;
        end
        ajbarcell{rr,qq} = tempajbarcell;
    end
end
etacell = cell(size(neffcellmmf));
for rr = 1:length(neffcellmmf)
    for qq = 1:length(neffcellmmf)
        if isempty(neffcellmmf{rr,qq})
            break
        end
        tempbetacell = betacell{rr,qq};
        tempajbarcell = ajbarcell{rr,qq};
        tempeta = zeros(1,length(lambdaarr));
        for mm = 1:length(lambdaarr)
            temparr = cell2mat(tempbetacell(1,mm));
            temparr1 = cell2mat(tempajbarcell(1,mm));
            for nn = 1:length(temparr)
                beth = temparr(nn);
                ah = temparr1(nn);
                for pp = 1:length(temparr)
                    tempeta(1,mm) = tempeta(1,mm) + (((temparr1(pp))^2)*(ah^2)*exp((1j)*(temparr(pp) - beth)*zmmf));
                end
            end
        end
        etacell{rr,qq} = tempeta;
    end
end

for nn = 1:length(neffcellmmf)
    figure(nn)
    hold on
    ylabel('Coupling Efficiency (dB)')
    xlabel('Wavelength (micron)')
    for mm = 1:length(neffcellmmf)
        if isempty(neffcellmmf{mm,nn})
            break
        end
        plot(lambdaarr,10.*log(cell2mat(etacell(mm,nn)))./log(10))
        xlim([1.5 1.6])
    end
end

% pklambda = zeros(size(neffcellmmf));
% for nn = 1:length(neffcellmmf)
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{mm,nn})
%             break
%         end
%         tempe = cell2mat(etacell(mm,nn));
%         pklambda(mm,nn) = lambdaarr(find(tempe == max(tempe)));
%     end
% end
% nummode = zeros(szneffcellmmf(1),szneffcellmmf(2));
% for aa = 1:szneffcellmmf(1)
%     for bb = 1:szneffcellmmf(2)
%         nummode(aa,bb) = length(cell2mat(neffcellmmf(aa,bb)));
%     end
% end
% plot(lambdaarr,nummode(1,:),'-x')
% hold on
% plot(lambdaarr,nummode(2,:),'-x')
% plot(lambdaarr,nummode(3,:),'-x')
% plot(lambdaarr,nummode(4,:),'-x')
% plot(lambdaarr,nummode(5,:),'-x')
% plot(lambdaarr,nummode(6,:),'-x')
% xlim([1.5 1.6])

% plot(0:0.1:0.5,[lambdaarr(find(etaarr(1,:) == max(etaarr(1,:)))) lambdaarr(find(etaarr(2,:) == max(etaarr(2,:)))) lambdaarr(find(etaarr(3,:) == max(etaarr(3,:)))) lambdaarr(find(etaarr(4,:) == max(etaarr(4,:)))) lambdaarr(find(etaarr(5,:) == max(etaarr(5,:)))) lambdaarr(find(etaarr(6,:) == max(etaarr(6,:))))],'-x')
% plot(0:0.1:0.5,[10*log(max(etaarr(1,:)))/log(10) 10*log(max(etaarr(2,:)))/log(10) 10*log(max(etaarr(3,:)))/log(10) 10*log(max(etaarr(4,:)))/log(10) 10*log(max(etaarr(5,:)))/log(10) 10*log(max(etaarr(6,:)))/log(10)],'-x')
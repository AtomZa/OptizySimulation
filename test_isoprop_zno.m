clear all
close all
clc
%-----------------------------------------------------------------------------%
global ncoresmf ncladdingsmf asmf ncoremmf ncladdingmmf ammf zmmf epscoresmf epscladdingsmf epscoremmf epscladdingmmf k0 neffarrsmf lambda nmedium aair
zno_data = load('neff_zno_data.mat');
isoprop_data = load('neff_isoprop_data_10%.mat');
volfact = zno_data.volfact;
lambda = 0.633; %wavelength
% lambda = 0.633;
ncoresmf = 1.4504; %core refractive index of smf
ncladdingsmf = 1.447; %cladding refractive index of smf
asmf = 4.1; %radius of smf
ncoremmf = 1.443; %refractive index of no core mmf
ammf = 62.5; %radius of no core mmf
zmmf = 57713; %length of mmf segment
epscoresmf = ncoresmf^2; %core pzermittivity of smf
epscladdingsmf = ncladdingsmf^2; %cladding permittivity of smf
epscoremmf = ncoremmf^2; %permittivity of no core mmf
k0 = 2*pi/lambda; %wave number
aair = 87.5; %radius of medium

zarr = 50000:10:60000; %array of propagation length

%change zarr if the code take too long to run

% nsolution = isoprop_data.neffir(:,find(isoprop_data.lambdair == lambda));
nsolution = isoprop_data.neff633;
% nzno = zno_data.nir(find(zno_data.lambdair == lambda));
nzno = zno_data.n633;
nmedium = zeros(length(nsolution),1);
for nn = 1:length(nsolution)
    for mm = 1:length(volfact)
        nmedium(nn,mm) = (volfact(mm)*nzno) + ((1-volfact(mm))*nsolution(nn));
    end
end
ncladdingmmf = nmedium; %refractive index of medium surrounding no core mmf
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
% plot(rarrall,eall)
%----------------------finding mode in mmf--------------------------------%
neffcellmmf = cell(length(nsolution),length(volfact));
for qq = 1:length(nsolution)
    for nn = 1:length(volfact)
        if ncladdingmmf(qq,nn) >= ncoremmf
            break
        end
        ifinal = ((ncoremmf - ncladdingmmf(qq,nn))/increment);
        epscladdingmmf = ncladdingmmf(qq,nn)^2; %permittivity of medium surrounding no core mmf
        tempmmf = zeros(1,round(ifinal)); %empty array of effective indices of mmf to be assigned in the next loop
        for ii = 0:ifinal
            try 
                chk = fzero(@characeqnmmf,[(ncladdingmmf(qq,nn) + (ii*increment)),(ncladdingmmf(qq,nn) + ((ii+1)*increment))]);
            catch
                tempmmf(ii+1) = 0;
             continue
            end
                tempmmf(ii+1) = fzero(@characeqnmmf,[(ncladdingmmf(qq,nn) + (ii*increment)),(ncladdingmmf(qq,nn) + ((ii+1)*increment))]);
        end
        tempmmf = tempmmf(tempmmf ~= 0);
        for p = 1:(length(tempmmf))/2
            temp = tempmmf(p);
            tempmmf(p) = tempmmf(length(tempmmf) - (p - 1));
            tempmmf(length(tempmmf) - (p - 1)) = temp;
        end
        neffcellmmf{qq,nn} = tempmmf;
    end
end
% neffcellmmf = neffcellmmf(~cellfun('isempty',neffcellmmf));
%------------------------finding u and w for each mode in mmf----------------------------------%
ucellmmf = cell(size(neffcellmmf)); %u of each mode in mmf
wcellmmf = cell(size(neffcellmmf)); %w of each mode in mmf
ccellmmf = cell(size(neffcellmmf)); %scaling coefficient
for qq = 1:length(neffcellmmf)
    for j = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,j})
            break
        end
        tempneffmmf = cell2mat(neffcellmmf(qq,j));
        tempu = zeros(size(tempneffmmf));
        tempw = zeros(size(tempneffmmf));
        tempc = zeros(size(tempneffmmf));
        epscladdingmmf = ncladdingmmf(j)^2;
        for jj = 1:length(tempneffmmf)
            tempu(jj) = findummf(tempneffmmf(jj));
            tempw(jj) = findwmmf(tempneffmmf(jj));
            tempc(jj) = (besselj(0,tempu(jj)))/(besselk(0,tempw(jj)));
        end
        ucellmmf{qq,j} = tempu;
        wcellmmf{qq,j} = tempw;
        ccellmmf{qq,j} = tempc;
    end
end
%------------------------finding field and its intensity for each mode in mmf----------------------------------%

ecellmmf = cell(size(neffcellmmf));
icellmmf = cell(size(neffcellmmf)); %intensity of each mode in mmf
for qq = 1:length(neffcellmmf)
    for nn = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,nn})
            break
        end
    	tempneffmmf = cell2mat(neffcellmmf(qq,nn));
        tempu = cell2mat(ucellmmf(qq,nn));
        tempw = cell2mat(wcellmmf(qq,nn));
        tempc = cell2mat(ccellmmf(qq,nn));
        earrallmmf = zeros(length(tempneffmmf),length(rarrall)); %field of each mode in mmf
        for k = 1:length(tempneffmmf)
            for m = 1:length(rarrall)
                if (rarrall(m) > ammf)
                    earrallmmf(k,m) = tempc(k)*(besselk(0,(tempw(k)*rarrall(m)/ammf)));
                elseif (rarrall(m) < -ammf)
                    earrallmmf(k,m) = tempc(k)*(besselk(0,(tempw(k)*rarrall(length(rarrall) - (m - 1))/ammf)));
                else
                    earrallmmf(k,m) = besselj(0,(tempu(k)*rarrall(m)/ammf));
                end
            end
        end
        ecellmmf{qq,nn} = earrallmmf;
        icellmmf{qq,nn} = (abs(earrallmmf)).^2;
    end
end
%-----------------------finding variable pj and aj for each mode in mmf-----------------------------------------%
pjcell = cell(size(neffcellmmf));
ajcell = cell(size(neffcellmmf));
for qq = 1:length(neffcellmmf)
    for mm = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,mm})
            break
        end
        tempneffmmf = cell2mat(neffcellmmf(qq,mm));
        tempu = cell2mat(ucellmmf(qq,mm));
        tempw = cell2mat(wcellmmf(qq,mm));
        tempc = cell2mat(ccellmmf(qq,mm));
        temppj = zeros(size(tempneffmmf)); %empty array of parameter pj 
        for nn = 1:length(tempneffmmf)
            funpj1 = @(r) ((abs(besselj(0,(tempu(nn).*r/ammf)))).^2).*r; %function of pj where r is between 0 and ammf
            funpj2 = @(r) ((abs(tempc(nn).*besselk(0,(tempw(nn).*r/ammf)))).^2).*r; %function of pj where r is greater than ammf
            temppj(nn) = ((integral(funpj1,0,ammf)) + (integral(funpj2,ammf,Inf))).*2.*pi;
        end
        pjcell{qq,mm} = temppj;
    end
end
for qq = 1:length(neffcellmmf)
    for mm = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,mm})
            break
        end
        tempneffmmf = cell2mat(neffcellmmf(qq,mm));
        tempu = cell2mat(ucellmmf(qq,mm));
        tempw = cell2mat(wcellmmf(qq,mm));
        tempc = cell2mat(ccellmmf(qq,mm));
        temppj = cell2mat(pjcell(qq,mm));
        tempaj = zeros(size(tempneffmmf)); %empty array of parameter aj
        integ1 = zeros(size(tempneffmmf)); %empty arrays for storing integrals in the next loop
        integ2 = zeros(size(tempneffmmf));
        integ3 = zeros(size(tempneffmmf));
        for nn = 1:length(tempneffmmf)
%     funaj1 = @(r) (cross((besselj(0,(usmf.*r./asmf))),conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;
%     funaj2 = @(r) (cross((besselk(0,(usmf.*r./asmf))),conj(besselj(0,(uarrmmf(nn).*r/ammf))))).*r;
%     funaj3 = @(r) (cross(eair2(length(eair2)),conj(besselk(0,(uarrmmf(nn).*r/ammf))))).*r;
% error with cross product so scalar multiplication is applied instead %
            funaj1 = @(r) ((besselj(0,(usmf.*r./asmf))).*(conj(besselj(0,(tempu(nn).*r/ammf))))).*r; %function where r is between 0 and asmf
            funaj2 = @(r) ((c.*besselk(0,(wsmf.*r./asmf))).*(conj(besselj(0,(tempu(nn).*r/ammf))))).*r; %function where r is between asmf and ammf
            funaj3 = @(r) ((eair2(length(eair2))).*(conj(tempc(nn).*besselk(0,(tempw(nn).*r/ammf))))).*r; %function where r is greater than ammf
            integ1(nn) = (integral(funaj1,0,asmf)); %integrate from 0 to asmf
            integ2(nn) = (integral(funaj2,asmf,ammf)); %integrate from asmf to ammf
            integ3(nn) = (integral(funaj3,ammf,Inf)); %integrate from ammf to inf
            tempaj(nn) = (integ1(nn) + integ2(nn) + integ3(nn)).*2.*pi./temppj(nn);
        end
        ajcell{qq,mm} = tempaj;
    end
end
%-----------------------finding efield of mmf with interference---------------------------------------%
betcell = cell(size(neffcellmmf));
for qq = 1:length(neffcellmmf)
    for nn = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,nn})
            break
        end
        tempneffmmf = cell2mat(neffcellmmf(qq,nn));
        betcell{qq,nn} = tempneffmmf.*k0; %array of propagation constants of each mode
    end
end
% eintercell = cell(size(neffcellmmf));
% iintercell = cell(size(neffcellmmf)); %array of field intensity of multimode interference
% einter0cell = cell(size(neffcellmmf));
% iinter0cell = cell(size(neffcellmmf));  %array of field intensity of multimode interference at z = 0
% for ss = 1:length(neffcellmmf)
%     for qq = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{ss,qq})
%             break
%         end
%         tempneffmmf = cell2mat(neffcellmmf(ss,qq));
%         tempbet = cell2mat(betcell(ss,qq));
%         tempaj = cell2mat(ajcell(ss,qq));
%         tempemmf = cell2mat(ecellmmf(ss,qq));
%         einter = zeros(length(zarr),length(rarrall)); %array of field of multimode interference
%         einter0 = zeros(1,length(rarrall)); %array of field of multimode interference at z = 0
%         for nn = 1:length(zarr)
%             for mm = 1:length(rarrall)
%                 for pp = 1:length(tempneffmmf)
%                     einter(nn,mm) = einter(nn,mm) + (tempaj(pp)*tempemmf(pp,mm)*exp(1j*tempbet(pp)*zarr(nn)));
%                 end
%             end
%         end
%         for mm = 1:length(rarrall)
%             for pp = 1:length(tempneffmmf)
%                 einter0(mm) = einter0(mm) + (tempaj(pp)*tempemmf(pp,mm)*exp(1j*tempbet(pp)*0));
%             end
%         end
%         eintercell{ss,qq} = einter;
%         iintercell{ss,qq} = (abs(einter)).^2;
%         einter0cell{ss,qq} = einter0;
%         iinter0cell{ss,qq} = (abs(einter0)).^2;
%     end
% end
% contourf(rarrall,zarr,iinter,'LineStyle','none','LevelList',0:0.01:1)
% hold on
% colormap(jet)

%-----------------------finding coupling efficiency of mmf with interference----------------------------------%
funescore = @(r) ((abs(besselj(0,(usmf.*r./asmf)))).^2).*r;
funesclad = @(r) ((abs(c*besselk(0,(wsmf.*r./asmf)))).^2).*r;
funesair = @(r) ((ecladdingsmf2(length(ecladdingsmf2))).^2).*r;
% ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + integral(funesair,ammf,Inf)).*(2*pi);
ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf) + 0).*(2*pi); %assume third integral is zero
ajbarcell = cell(size(neffcellmmf));
for qq = 1:length(neffcellmmf)
    for nn = 1:length(neffcellmmf)
        if isempty(neffcellmmf{qq,nn})
            break
        end
        tempaj = cell2mat(ajcell(qq,nn));
        temppj = cell2mat(pjcell(qq,nn));
        ajbarcell{qq,nn} = tempaj.*(sqrt(temppj./ps)); %aj bar parameters of each mode
    end
end

etacell = cell(size(neffcellmmf));
for ss = 1:length(neffcellmmf)
    for qq = 1:length(neffcellmmf)
        if isempty(neffcellmmf{ss,qq})
            break
        end
        tempneffmmf = cell2mat(neffcellmmf(ss,qq));
        tempbet = cell2mat(betcell(ss,qq));
        tempajbar = cell2mat(ajbarcell(ss,qq));
        etaarr = zeros(1,length(zarr));
        for nn = 1:length(zarr)
            for mm = 1:length(tempneffmmf)
                beth = tempbet(mm);
                ahb = tempajbar(mm);
                for pp = 1:length(tempneffmmf)
                    ajb = tempajbar(pp);
                    betj = tempbet(pp);
                    etaarr(1,nn) = etaarr(1,nn) + ((ajb^2)*((conj(ahb))^2)*exp(1j*(betj - beth)*zarr(nn)));
                end
            end
        end
        etacell{ss,qq} = etaarr;
    end
end
% eta0cell = cell(size(neffcellmmf));
% for ss = 1:length(neffcellmmf)
%     for qq = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{ss,qq})
%             break
%         end
%         tempneffmmf = cell2mat(neffcellmmf(ss,qq));
%         tempbet = cell2mat(betcell(ss,qq));
%         tempajbar = cell2mat(ajbarcell(ss,qq));
%         eta0 = 0;
%         for mm = 1:length(tempneffmmf)
%             beth = tempbet(mm);
%             ahb = tempajbar(mm);
%             for pp = 1:length(tempneffmmf)
%                 ajb = tempajbar(pp);
%                 betj = tempbet(pp);
%                 eta0 = eta0 + ((ajb^2)*((conj(ahb))^2)*exp(1j*(betj - beth)*0));
%             end
%         end
%         eta0cell{ss,qq} = eta0;
%     end
% end
% 
% coupcell = cell(size(neffcellmmf));
% zimcell = cell(size(neffcellmmf));
% for nn = 1:length(neffcellmmf)
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{nn,mm})
%             break
%         end
%         tempeta = cell2mat(etacell(nn,mm));
%         maxe = max(tempeta);
%         coupcell{nn,mm} = real(10*log(maxe)/log(10));
%         zimcell{nn,mm} = zarr(find(tempeta == maxe));
%     end
% end
% 
% for nn = 1:length(neffcellmmf)
%     tempeta = zeros(1,length(neffcellmmf));
%     tempneff = zeros(1,length(neffcellmmf));
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{nn,mm})
%             break
%         end
%         tempeta(1,mm) = cell2mat(coupcell(nn,mm));
%         tempneff(1,mm) = nmedium(nn,mm);
%     end
%     tempeta = tempeta(tempeta ~= 0);
%     tempneff = tempneff(tempneff ~= 0);
%     tempvolfact = volfact(1,1:(mm-1));
%     figure(nn)
%     yyaxis('left')
%     plot(tempvolfact,tempeta,'-x')
%     ylabel('Coupling Efficiency (dB)')
%     yyaxis('right')
%     plot(tempvolfact,tempneff,'-x')
%     ylabel('Effective Index')
%     xlabel('ZnO Volume Fraction')
%     xlim([min(tempvolfact) max(tempvolfact)]) 
% end    
% 
% for qq = 1:length(neffcellmmf)
%     tempneff = zeros(1,length(neffcellmmf));
%     tempzim = zeros(1,length(neffcellmmf));
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{qq,mm})
%             break
%         end
%         tempneff(1,mm) = nmedium(qq,mm);
%         tempzim(1,mm) = cell2mat(zimcell(qq,mm));
%     end
%     tempneff = tempneff(tempneff ~= 0);
%     tempvolfact = volfact(1,1:(mm-1));
%     tempzim = tempzim(tempzim ~= 0);
%     figure(qq + nn)
%     yyaxis('left')
%     plot(tempvolfact,tempzim,'-x')
%     ylabel('Reimaging Distance (micron)')
%     yyaxis('right')
%     plot(tempvolfact,tempneff,'-x')
%     ylabel('Effective Index')
%     xlabel('ZnO Volume Fraction')
%     xlim([min(tempvolfact) max(tempvolfact)]) 
% end    
% 
% for ss = 1:length(neffcellmmf)
%     tempneff = zeros(1,length(neffcellmmf));
%     tempzim = zeros(1,length(neffcellmmf));
%     counter = 0;
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{mm,ss})
%             break
%         end
%         tempneff(1,mm) = nmedium(mm,ss);
%         tempzim(1,mm) = cell2mat(zimcell(mm,ss));
%         counter = counter + 1;
%     end
%     tempneff = tempneff(tempneff ~= 0);
%     tempvolfact = volfact(1,1:counter);
%     tempzim = tempzim(tempzim ~= 0);
%     figure(ss + qq + nn)
%     yyaxis('left')
%     plot(tempvolfact,tempzim,'-x')
%     ylabel('Reimaging Distance (micron)')
%     yyaxis('right')
%     plot(tempvolfact,tempneff,'-x')
%     ylabel('Effective Index')
%     xlabel('Isopropanol Volume Fraction')
%     xlim([min(tempvolfact) max(tempvolfact)]) 
% end 
% 
% for rr = 1:length(neffcellmmf)
%     tempneff = zeros(1,length(neffcellmmf));
%     tempeta = zeros(1,length(neffcellmmf));
%     counter = 0;
%     for mm = 1:length(neffcellmmf)
%         if isempty(neffcellmmf{mm,rr})
%             break
%         end
%         tempneff(1,mm) = nmedium(mm,rr);
%         tempeta(1,mm) = cell2mat(coupcell(mm,rr));
%         counter = counter + 1;
%     end
%     tempneff = tempneff(tempneff ~= 0);
%     tempvolfact = volfact(1,1:counter);
%     tempeta = tempeta(tempeta ~= 0);
%     figure(ss + qq + nn + rr)
%     yyaxis('left')
%     plot(tempvolfact,tempeta,'-x')
%     ylabel('Coupling Efficiency (dB)')
%     yyaxis('right')
%     plot(tempvolfact,tempneff,'-x')
%     ylabel('Effective Index')
%     xlabel('Isopropanol Volume Fraction')
%     xlim([min(tempvolfact) max(tempvolfact)]) 
% end 

for mm = 1:length(neffcellmmf)
    figure(mm)
    for nn = 1:length(neffcellmmf)
        if isempty(neffcellmmf{nn,mm})
            break
        end
        hold on
        plot(zarr,10.*log(cell2mat(etacell(nn,mm)))./log(10));
        xlabel('Porpagation Distance (micron)')
        ylabel('Coupling Efficiency (dB)')
    end
end
% % plot(zarr,10.*log(etaarr(1,:))./log(10))
% 
% % zim = [zarr(find(etaarr(1,:) == max(etaarr(1,:)))) zarr(find(etaarr(2,:) == max(etaarr(2,:)))) zarr(find(etaarr(3,:) == max(etaarr(3,:)))) zarr(find(etaarr(4,:) == max(etaarr(4,:)))) zarr(find(etaarr(5,:) == max(etaarr(5,:)))) zarr(find(etaarr(6,:) == max(etaarr(6,:))))];
% % coup = [10*log(max(etaarr(1,:)))/log(10) 10*log(max(etaarr(2,:)))/log(10) 10*log(max(etaarr(3,:)))/log(10) 10*log(max(etaarr(4,:)))/log(10) 10*log(max(etaarr(5,:)))/log(10) 10*log(max(etaarr(6,:)))/log(10)];
% % yyaxis('left')
% % plot(0:0.1:0.5,zim,'-x')
% % hold on
% % yyaxis('right')
% % plot(0:0.1:0.5,coup,'-x')
% 
% % yyaxis('left')
% % plot(0:0.1:0.5,coup,'-x')
% % yyaxis('right')
% % plot(0:0.1:0.5,nmedium(1:6),'-x')
% 
% % yyaxis('left')
% % plot(0:0.1:0.5,[84 76 68 59 46 27],'-x')
% % hold on
% % yyaxis('right')
% % plot(0:0.1:0.5,nmedium(1:6,1),'-x')
% 
% % plot(volfact,nmedium,'-x')
% % y = ncoremmf.*(ones(size(volfact)));
% % hold on
% % plot(volfact,y,'--')
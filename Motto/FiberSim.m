clc
clearvars -except NSiO2 IPAindex ZnOindex
%% importing variable
load("matlab_simoutput\IPA(liquid)_(0.3-1.0).mat")
load("matlab_simoutput\SiO2(FusedSilica)_(0.21-1.5).mat")
load("matlab_simoutput\ZnO(ThinFilm)_(0.3-1.0).mat")
%% medium parameter
wavelength      = (0.3:0.01:1.0);    %wavelength array for testing (important)**
NSiO2           = NSiO2(1:791);      %wavelength dependent refractive index 
IPAindex        = IPAindex(1:10:701);%'
ZnOindex        = ZnOindex(1:10:701);%'
aair            = 72.5;              %radius of medium
aNP             = 0.5;               % Nanoparticle layers tickness
%% fiber parameter
asmf            = 4.1;               %radius of smf
ncoresmf        = 1.4504;            %refractive index of smf core
ncladdingsmf    = 1.447;             %refractive index  of smf cladding
ammf            = 62.5;              %radius of core mmf
ncoremmf        = NSiO2;             %refractive index of core mmf
zarr            = (0:50:150000);     %Propagation distant array     (important)**
NPsAlpha        = 0.0;               %NPs attenuation coefficient
mmfmode         = 80;                %Number of modes used
%% in-Loop parameter
rarrcosmf   = 0   :0.01:asmf;
rarrclsmf   = asmf:0.01:ammf; 
rarrair     = ammf:0.05:aair-0.05;
rarrcommf   = 0   :0.05:ammf-0.05;
rarrclmmf   = rarrair;
rarrmmf     = [rarrcommf rarrclmmf];
Radiusarr   = [flip(rarrmmf.*-1) rarrmmf];
dx          = 0.0001;
ifinal      = ((ncoresmf - ncladdingsmf)/dx);
usmf        = zeros(1,length(wavelength));
wsmf        = zeros(1,length(wavelength));
csmf        = zeros(1,length(wavelength));
eair        = zeros(1,length(wavelength));
eta         = zeros(length(zarr),length(wavelength));
% efil        = zeros(length(zarr),length(rarrmmf),length(wavelength));
nmedium = zeros(100,length(wavelength));
for i = 1:100
    nmedium = 
end

%% numerical sim section

for WL = 1:length(wavelength)   %SMF
    wavenum = 2*pi/wavelength(WL);
    characeqnsmf = @(neff) characeqn(neff,ncoresmf,ncladdingsmf,wavenum,asmf);
    for i = 0:ifinal
      try 
          chk = fzero(characeqnsmf,[(ncoresmf-((i+1)*dx)),ncoresmf-(i*dx)]);
          neffsmf = chk;
          break;
      catch 
          continue
      end
    end

    usmf(WL) =  wavenum*asmf*sqrt(ncoresmf^2 - neffsmf^2);   
    wsmf(WL) =  wavenum*asmf*sqrt(neffsmf^2 - ncladdingsmf^2);

    csmf(WL) = (besselj(0,usmf(WL)))/(besselk(0,wsmf(WL)));          %scaling coefficient
    eclsmf   = csmf(WL).*(besselk(0,(wsmf(WL).*rarrclsmf./asmf)));   %field where asmf < r < ammf
    ecosmf   = (besselj(0,(usmf(WL).*rarrcosmf./asmf)));
    eair(WL)   = csmf(WL).*(besselk(0,(wsmf(WL).*rarrair(end)./asmf)));  
end

parfor WL = 1:length(wavelength)% MMF
    wavenum = 2*pi/wavelength(WL);
    for NMED = 1:length(nmedium)  %MEDIUM LOOP
        ncladdingmmf = nmedium(NMED,WL);
        wavenum = 2*pi/wavelength(WL);
        neffarrmmf = zeros(1,mmfmode);
        dx = 5*(10^-6);
        mode = 0;
        pointer = ncoremmf(WL);
        ifinal = ((ncoremmf(WL) - ncladdingmmf)/dx);
        characeqnmmf = @(neff) characeqn(neff,ncoremmf(WL),ncladdingmmf,wavenum,ammf);
        for i = 0:ifinal
            try 
                chk = fzero(characeqnmmf,[pointer-dx,pointer]);
                mode = mode+1;
                neffarrmmf(mode) = chk;
            catch
                pointer = pointer-dx;
                if(pointer-dx < ncladdingmmf)
                    break;
                end
                continue
            end
            pointer = pointer-dx;
            if(mode>= mmfmode)
               break;
            end
            if(mode > 1)
            dx = abs(neffarrmmf(mode-1)-neffarrmmf(mode));
            end
        end
        neffarrmmf = neffarrmmf(neffarrmmf ~= 0);

%------------------------finding u and w for each mode in mmf----------------------------------%
    uarrmmf = zeros(1,mmfmode); %u of each mode in mmf
    warrmmf = zeros(1,mmfmode); %w of each mode in mmf
    carr = zeros(1,mmfmode); %scaling coefficient
    for j = 1:length(neffarrmmf)
        uarrmmf(j) =  wavenum.*ammf.*sqrt(ncoremmf(WL).^2 - neffarrmmf(j).^2);   
        warrmmf(j) =  wavenum.*ammf.*sqrt(neffarrmmf(j).^2 - ncladdingmmf.^2);
        carr(j) = (besselj(0,uarrmmf(j)))/(besselk(0,warrmmf(j),1));
    end
%-----------------------finding variable pj and aj for each mode in mmf-------------------------%
    pj = zeros(1,mmfmode); %empty array of parameter pj 
    aj = zeros(1,mmfmode); %empty array of parameter aj

    for ii = 1:length(neffarrmmf)
        funpj1 = @(r) ((abs(besselj(0,(uarrmmf(ii).*r/ammf)))).^2).*r;           %function of pj where r is between 0 and ammf
        funpj2 = @(r) ((abs(carr(ii).*besselk(0,(warrmmf(ii).*r/ammf)))).^2).*r; %function of pj where r is greater than ammf
        pj(ii) = ((integral(funpj1,0,ammf)) + (integral(funpj2,ammf,ammf*3))).*2.*pi;

        funaj1 = @(r) ((besselj(0,(usmf(WL).*r./asmf))).*(conj(besselj(0,(uarrmmf(ii).*r/ammf))))).*r;      %function where r is between 0 and asmf
        funaj2 = @(r) ((csmf(WL).*besselk(0,(wsmf(WL).*r./asmf))).*(conj(besselj(0,(uarrmmf(ii).*r/ammf))))).*r; %function where r is between asmf and ammf
        funaj3 = @(r) ((eair(WL)).*(conj(carr(ii).*besselk(0,(warrmmf(ii).*r/ammf))))).*r;                    %function where r is greater than ammf

        aj(ii) = (integral(funaj1,0,asmf) + integral(funaj2,asmf,ammf) + integral(funaj3,ammf,ammf*3)).*2.*pi./pj(ii);
    end
%--------------------------------finding Zeta(Total vs.clading power ratio) -------------------------------------------------%
      zetammf = zeros(1,mmfmode);
%     for iii = 1:length(neffarrmmf)
%         Jfunc   = @(r) besselj(0,uarrmmf(iii).*r./ammf);
%         Kfunc   = @(r) besselk(0,warrmmf(iii).*r./ammf,1);
%         cladfun = @(r) (carr(iii).*Kfunc(r).*exp((1-r./ammf).*warrmmf(iii))).^2;        
%         corefun = @(r) Jfunc(r).^2; 
%     
%         zetammf(iii) = integral(cladfun,ammf,ammf+aNP)./(integral(cladfun,ammf,aair) + integral(corefun,0,ammf));
%     end

 %-----------------------finding efield of mmf with interference---------------------------------------%

% earrallmmf = zeros(length(neffarrmmf),length(rarrmmf)); %field of each mode in mmf
% 
% for mode2=1:length(neffarrmmf)
%     earrallmmf(mode2,:) = [besselj(0,rarrcommf.*uarrmmf(mode2)./ammf) carr(mode2).*besselk(0,rarrclmmf.*warrmmf(mode2)./ammf)];
% end    

betammf = (neffarrmmf.*wavenum); %array of propagation constants of each mode% 
% zarr    = gpuArray(zarr);    %Gpu
% rarrmmf = gpuArray(rarrmmf); %Gpu
% nn      = int16(1:length(zarr));
% mm      = int16(1:length(rarrmmf));
% 
% p1 = (earrallmmf(:,mm)).*aj.';
% p2 = exp(((1j.*(betammf.'))-((zetammf.').*Alpha./2) ).*(zarr(nn))); 
% 
% einter = (p2.')*p1;
% iinter = gather(abs(einter)).^2; %array of field intensity of multimode interference
% einter = gather(einter.');
% 
% efil(:,:,WL) = iinter;
% rescaledMatrix = uint8(rescale(iinter, 0, 255));
% doubleRe = flip([flip(iinter,2) iinter],2);
% doubleRe = flip([flip(efil(:,:,8),2) efil(:,:,8)],2);
% imagesc(Radiusarr,zarr,doubleRe)
% xlabel('radius(um)')
% ylabel('Propagation distant(um)')
% title('Eletric field Wavelength = 0.763 um')
% colorbar
 %-----------------------finding coupling efficiency of mmf with interference----------------------------------%

    funescore = @(r) ((abs(besselj(0,(usmf(WL).*r./asmf)))).^2).*r;
    funesclad = @(r) ((abs(csmf(WL)*besselk(0,(wsmf(WL).*r./asmf)))).^2).*r;
    ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf)).*(2*pi); %assume third integral is zero
    ajbar = aj.*(sqrt(pj./ps));
    
%     ======================== (ZARR,WAVELENGT)============================
%     nn = 1:length(zarr);
% 
%     P1= reshape((((betammf - betammf.')*1j)-((zetammf - zetammf.').*(Alpha/2))).',[],1).*zarr(nn);
%     P2= reshape(log((ajbar.^2).*((conj(ajbar)).^2).'),1,[]);
%     P3 = sum(exp(P2+P1.'),2);
% 
%     eta(:,WL) =  abs(reshape(sum(P3,2),[],1)).';
%     =========================(WAVELENGTH,NMED)=========================== 
       P1 = (betammf - betammf.')*1j.*(103000);   
       P2 = log((ajbar.^2).*((conj(ajbar)).^2).');
       P3 = sum(exp(P2+P1.'));
       eta(NMED,WL) =  abs(reshape(sum(P3,2),[],1)).';
%     =====================================================================
    end
end
clearvars -except efil wavelength zarr rarrmmf Radiusarr eta
%% function 
function f = characeqn(neff,ncore,nclad,k,a)
u = k*a*sqrt(ncore.^2 - neff.^2);
w = k*a*sqrt(neff.^2 - nclad.^2);
f = w*besselj(0,u)*besselk(1,w,1)-u*besselk(0,w,1)*besselj(1,u);
end

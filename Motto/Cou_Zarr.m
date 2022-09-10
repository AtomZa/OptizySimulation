clc
clearvars -except
%% medium parameter
wavelength      = (0.3:0.0005:0.9);   % Wavelength array for testing  (important)**
aair            = 72.5;               % Radius of medium + mmf_core
aNP             = 0.5;                % NPs layers tickness
NPsAlpha        = 0;                  % NPs attenuation coefficient
fVOC            = 0;                  % VOC Volume fraction           (important)**
fZnO            = 0;                  % ZnO Volume fraction           (important)**
NSiO2           = SiO2_f(wavelength);   
IPAindex        = IPA_f(wavelength);
ACNindex        = ACN_f(wavelength);
EtOHindex       = EtOH_f(wavelength);
Metindex        = Met_f(wavelength);
ZnOindex        = ZnO_f(wavelength); 
%% fiber parameter
asmf            = 4.1;               % Radius of SMF
ncoresmf        = 1.4504;            % Refractive index of SMF core
ncladdingsmf    = 1.447;             % Refractive index of SMF cladding
ammf            = 62.5;              % Radius of core MMF
ncoremmf        = NSiO2;             % Refractive index of core MMF
zarr            = (0:20:50000);     % Propagation distant array     (important)**
mmfmode         = 40;                % Number of radial modes used

% param = saveparam(ammf,asmf,0.0005,50,fIPA,fZnO);
% exportT(param,'Coupling_condition_A.csv')
%% in-Loop parameter
rarrcosmf   = 0:0.01:asmf;
rarrclsmf   = asmf:0.01:ammf; 
rarrair     = ammf:0.05:aair-0.05;
dx          = 0.0001;
ifinal      = ((ncoresmf - ncladdingsmf)/dx);
usmf        = zeros(1,length(wavelength));
wsmf        = zeros(1,length(wavelength));
csmf        = zeros(1,length(wavelength));
eair        = zeros(1,length(wavelength));
eta         = zeros(length(zarr),length(wavelength)); 
nmedium      = effmedium(IPAindex,ZnOindex,fVOC,fZnO);  % Create medium refractive index array
%% Simulation
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
        ncladdingmmf = nmedium(WL);
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
    uarrmmf = zeros(1,length(neffarrmmf)); %u of each mode in mmf
    warrmmf = zeros(1,length(neffarrmmf)); %w of each mode in mmf
    carr = zeros(1,length(neffarrmmf)); %scaling coefficient
    for j = 1:length(neffarrmmf)
        uarrmmf(j) =  wavenum.*ammf.*sqrt(ncoremmf(WL).^2 - neffarrmmf(j).^2);   
        warrmmf(j) =  wavenum.*ammf.*sqrt(neffarrmmf(j).^2 - ncladdingmmf.^2);
        carr(j) = (besselj(0,uarrmmf(j)))/(besselk(0,warrmmf(j),1));
    end
%-----------------------finding variable pj and aj for each mode in mmf-------------------------%
    pj = zeros(1,length(neffarrmmf)); %empty array of parameter pj 
    aj = zeros(1,length(neffarrmmf)); %empty array of parameter aj

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
%       zetammf = zeros(1,length(neffarrmmf));
%     for iii = 1:length(neffarrmmf)
%         Jfunc   = @(r) besselj(0,uarrmmf(iii).*r./ammf);
%         Kfunc   = @(r) besselk(0,warrmmf(iii).*r./ammf,1);
%         cladfun = @(r) (carr(iii).*Kfunc(r).*exp((1-r./ammf).*warrmmf(iii))).^2;        
%         corefun = @(r) Jfunc(r).^2; 
%     
%         zetammf(iii) = integral(cladfun,ammf,ammf+aNP)./(integral(cladfun,ammf,aair) + integral(corefun,0,ammf));
%     end
 %-----------------------finding coupling efficiency of mmf with interference----------------------------------%
    betammf = (neffarrmmf.*wavenum); %array of propagation constants of each mode% 
    funescore = @(r) ((abs(besselj(0,(usmf(WL).*r./asmf)))).^2).*r;
    funesclad = @(r) ((abs(csmf(WL)*besselk(0,(wsmf(WL).*r./asmf)))).^2).*r;
    ps = (integral(funescore,0,asmf) + integral(funesclad,asmf,ammf)).*(2*pi); %assume third integral is zero
    ajbar = aj.*(sqrt(pj./ps));
    
%     ======================== (ZARR,WAVELENGT)============================
  
%  P1= reshape((((betammf - betammf.')*1j)-((ones(length(neffarrmmf)).*zetammf).*(NPsAlpha/2))).',[],1).*zarr;
   P1= reshape((((betammf - betammf.')*1j)).',[],1).*zarr;
   P2= reshape(log((ajbar.^2).*((conj(ajbar)).^2).'),1,[]);
   P3 = sum(exp(P2+P1.'),2);

   eta(:,WL) =  abs(P3);
end

% clearvars -except wavelength zarr eta 
clearvars -except efil wavelength zarr rarrmmf Radiusarr eta betammf zetammf
imagesc(zarr(100:end),wavelength,eta(100:end,:).')
set(gca,'YDir','normal')
xlabel('Propagation Distant (um)')
ylabel('Wavelength (um)')
grid
grid minor

%% function
function f = characeqn(neff,ncore,nclad,k,a)
u = k*a*sqrt(ncore.^2 - neff.^2);
w = k*a*sqrt(neff.^2 - nclad.^2);
f = w*besselj(0,u)*besselk(1,w,1)-u*besselk(0,w,1)*besselj(1,u);
end

function f = effmedium(nVoc,nCoating,fVoc,fCoating)
nsol = (fVoc/100).*nVoc(:)+((1-(fVoc/100))*1.0003);  %change 1 to 0.06374 for realism
f    = ((fCoating/100).*nCoating(:))+((1-(fCoating/100)).*nsol);
end

function f = SiO2_f(lambda) %Malitson 1965: Fused silica; n 0.21-6.7 µm
f1 =  (0.6961663.*(lambda.^2))./((lambda.^2)-(0.0684043^2));
f2 =  (0.4079426.*(lambda.^2))./((lambda.^2)-(0.1162414^2));
f3 =  (0.8974794.*(lambda.^2))./((lambda.^2)-(9.896161^2));
f  = sqrt(f1+f2+f3+1);
end

function f = IPA_f(lambda) %Sani and Dell'Oro 2016: iso-propanol;Sellmier; n,k 0.185-2.8 µm
f1 =  (0.0107.*(lambda.^2))./((lambda.^2)- 8.88 );
f2 =  (0.8702.*(lambda.^2))./((lambda.^2)- 0.01036 );
f  = sqrt(f1+f2+1);
end

function f = ZnO_f(lambda)     %Y.Yang 2006:ZnO_ThinflimOnGlass; Cauchy; 0.4-0.8 µm 
f1 =  0.059./(lambda.^2);
f  = f1+1.781;
end

function f = ACN_f(lambda)     %J. Rheims, J Köser:Acetone_Abbe; --; 0.4765–0.83µm
f1 =  0.00306./(lambda.^2);
f2 =  0.00006./(lambda.^4);
f  = 1.34979+f1+f2;
end

function f = EtOH_f(lambda)    %Sani and Dell'Oro 2016: Ethanol;Sellmier; n,k 0.185-2.8 µm
f1 =  (0.0165.*(lambda.^2))./((lambda.^2)-(9.08));
f2 =  (0.8268.*(lambda.^2))./((lambda.^2)-(0.01039));
f  = sqrt(f1+f2+1);
end

function f = Met_f(lambda)     %Moutzouris et al. 2013; n 0.450- 1.55um
f1 =  1.745946239;
f2 =  (0.005362181.*(lambda.^2));
f3 =  (0.004656355.*(lambda.^-2));
f4 =  (0.00044714.*(lambda.^-4));
f5 =  (0.00001508.*(lambda.^-6));
f = sqrt(f1-f2+f3+f4-f5);
end

% function f = saveparam(ammf,asmf,res_wavelength,res_zarr,volmed,volcoat)
%  f = table(['mmf_radius';'smf_radius';'wavelength_resolution';'zarr_resolution'; ...
%             'Medium_VolumeFraction';'Coating_VolumeFraction'] ...
%             ,[ammf;asmf;res_wavelength;res_zarr;volmed;volcoat]);
% end

% function exportT(t,name)
%  writetable(t,name)
% end

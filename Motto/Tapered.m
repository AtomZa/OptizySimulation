clc
clearvars -except
%% medium parameter
wavelength      = 1;   %wavelength array for testing  (important)**
% aNP             = 0.5;               % NPs layers tickness
% NPsAlpha        = 0;                 % NPs attenuation coefficient
fIPA            = 0;                 % IPA Volume fraction           (important)**
fZnO            = 0;                 % ZnO Volume fraction           (important)**
NSiO2           = SiO2_f(wavelength);   
IPAindex        = IPA_f(wavelength);
ZnOindex        = ZnO_f(wavelength); 
nmedium      = effmedium(IPAindex,ZnOindex,fIPA,fZnO);  % Create medium refractive index array
%% fiber parameter
ncoresmf        = 1.450;             % Refractive index of smf core
ncladdingsmf    = 1.447;             % Refractive index of smf cladding
asmf            = 4.1;               % Radius of core smf (um)
ncoremmf        = NSiO2;             % Refractive index of core mmf
mmfmode         = 20;                % Number of modes used
dx              = 5*(10^-6);         % refractive index resolution

dZ              = 400;               % Taper resulotion (um)
waistL          = 10000;             % Taper waist length (um)
TransL          = 5000;              % Taper Transition length (um)
StartR          = 62.5;              % Taper input fiber radius (um)
WaistR          = 40;                % Taper waist radius (um)

Taper           = ( WaistR:(dZ*(abs(WaistR-StartR)/TransL)):StartR);
segment_R       = [ asmf flip(Taper) Taper asmf];

TotalPower      = zeros(1,length(segment_R));
a  = zeros(length(segment_R),mmfmode);     %empty array of parameter aj

%% Simulation

for WL = 1:length(wavelength) 
    k0 = 2*pi/wavelength(WL);
    ncladmmf = nmedium(WL);
    for seg = 1:length(segment_R) % Fiber segment
            if(seg == 1 || seg == length(segment_R))
                mmfmode = 1;
            else
                mmfmode = 20;
            end
            dx = 5*(10^-6); 
            ammfC = segment_R(seg);
            neff = zeros(1,mmfmode);
            pointer = ncoremmf(WL);
            ifinal = ((ncoremmf(WL) - ncladmmf)/dx);
            characeqnmmf = @(neff) characeqn(neff,ncoremmf(WL),ncladmmf,k0,ammfC);
            mode = 0;
            for i = 0:ifinal
                try 
                    chk = fzero(characeqnmmf,[pointer-dx,pointer]);
                    mode = mode+1;
                    neff(mode) = chk;
                catch
                    pointer = pointer-dx;
                    if(pointer-dx < ncladmmf)
                        break;
                    end
                    continue
                end
                pointer = pointer-dx;
                if(mode>= mmfmode)
                   break;
                end
                if(mode > 1)
                dx = abs(neff(mode-1)-neff(mode));
                end
            end
            neff = neff(neff ~= 0);
            b = k0*neff;
           %------------------------finding u and w of each mode ----------------------------------%
            u = zeros(1,length(neff)); 
            w = zeros(1,length(neff)); 
            c = zeros(1,length(neff)); % scaling coefficient
            for j = 1:length(neff)
                u(j) =  k0.*ammfC.*sqrt(ncoremmf(WL).^2 - neff(j).^2);   
                w(j) =  k0.*ammfC.*sqrt(neff(j).^2      - ncladmmf.^2);
                c(j) = (besselj(0,u(j)))/(besselk(0,w(j),1));
            end
%-----------------------finding variable pj and aj for each mode in mmf-------------------------%
        if(seg > 1 && seg ~= length(segment_R)) % smf for first and last
          ammfP = segment_R(seg-1);
          p  = zeros(1,length(neff));     %empty array of parameter pj 
%           a  = zeros(1,length(neff));     %empty array of parameter aj
          Rarr = linspace(0,ammfC*1.2,(u(end)*20*pi)*1.1);
            for j =  1:length(neff)  
                R_P = Efield_Discrete(up(j),wp(j),cp(j),ammfP,Rarr);
  
                if (j == 1) 
                    E_P = R_P;  
                else
                    if(ammfP==ammfC)
                        dL = waistL;
                    else
                        dL = dZ;
                    end
                    E_P = E_P+a(seg-1,j).*R_P.*exp(1i*dL*bp(j));
                    if(mod(j,2)==0 && seg==3) ,hold on,plot(Rarr,R_C);end
                end
            end
            
            for j =  1:length(neff)  
%                 hold on
                R_C     = Efield_Discrete(u(j),w(j),c(j),ammfC,Rarr);
                 if(mod(j,20)==0 && seg== 2),plot(Rarr,R_C);end
                p(j)    = Epower(u(j),w(j),c(j),ammfC);  
                a(seg,j)    = Expansion_aj_Discrete(E_P,R_C,p(j),ammfC,Rarr);
%                 TotalPower(seg) = TotalPower(seg)+p(j);
            end
            up = u;
            wp = w;
            bp = b;
            cp = c;
        elseif(seg == 1)
            up = u.*ones(1,20);
            wp = w.*ones(1,20);
            bp =   zeros(1,20);
            cp = c.*ones(1,20);
            a(seg,:) = zeros(1,20);
            a(seg,1) = 1;
        else
        end
    end
end
% plot(Rarr,R_C)
% hold on
% plot(Rarr,E_P)
% xline(ammfC)

% clearvars -except wavelength betammf

%% function //Appendix As
function f = Epower(u,w,c,R)
funpj1 = @(r) ((abs(besselj(0,(u.*r/R)))).^2).*r;            %function of pj where r is between 0 and ammf2
funpj2 = @(r) ((abs(c.*besselk(0,(w.*r/R)))).^2).*r;         %function of pj where r is greater than ammf2
f = ((integral(funpj1,0,R)) + (integral(funpj2,R,R*1.2))).*2.*pi;
end

function f = Efield_Discrete(u,w,c,R,Rarr)
funcore = @(r)   (r<R & r>=0).*besselj(0,(u.*r)/R);         
funclad = @(r)   (r>=R).*c.*besselk(0,(w.*r)/R);
index = find(0.1>abs(Rarr-R));

f = [funcore(Rarr(1:index)) funclad(Rarr(index+1:end))];
end

function f = Expansion_aj_Discrete(EP,EC,power,R,Rarr) 
 Y  = EP.*EC.*Rarr;
 dR = (R*1.2)/length(Rarr);
 f  = (dR*trapz(Y)*2*pi)/power;
end

function f = characeqn(neff,ncore,nclad,k,a)
u = k*a*sqrt(ncore.^2 - neff.^2);
w = k*a*sqrt(neff.^2 - nclad.^2);
f = w*besselj(0,u)*besselk(1,w,1)-u*besselk(0,w,1)*besselj(1,u);
end


function f = effmedium(nVoc,nCoating,fVoc,fCoating)
nsol = (1*fVoc/100).*nVoc(:)+((1-(1*fVoc/100))*1.0003);  %change 1 to 0.06374 for realism
f    = ((fCoating/100).*nCoating(:))+((1-(fCoating/100)).*nsol);
end

%% Refractive index (Dispersion relation) function //Appendix B
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

function f = ZnO_f(lambda) %Yang:ZnOthinflim_2006; Cauchy; Unknown range
f1 =  0.059./(lambda.^2);
f  = f1+1.781;
end

%% Un-use prototype function 
function f = Efield_Function(u,w,c,Core_radius)
funpj = @(r)   besselj(0,(u.*r./Core_radius)).^2.*r + c*besselk(0,(w.*r./Core_radius)).^2.*r;          %function of pj where r is between 0 and ammf
f = funpj;
end

function f = Expansion_aj(Ep,u,w,c,ammfP,ammf,power) % what it suppose to be

Faj1 = @(r) Ep(r).*(conj(besselj(0,(u.*r/ammf)))).*r;       %function where r is between 0 and ammf1
Faj2 = @(r) Ep(r).*(conj(besselj(0,(u.*r/ammf)))).*r;       %function where r is between ammf1 and ammf2
Faj3 = @(r) Ep(r).*(conj(c.*besselk(0,(w.*r/ammf)))).*r;   %function where r is greater than ammf2

f = (integral(Faj1,0,ammfP) + integral(Faj2,ammfP,ammf) + integral(Faj3,ammf,ammf*2)).*2.*pi./power;
end


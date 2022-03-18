function f = findusmf(neff)
global asmf k0 epscoresmf
f = k0.*asmf.*(sqrt(epscoresmf - (neff.^2)));
end

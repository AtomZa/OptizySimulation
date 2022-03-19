function f = findwsmf(neff)
global k0 asmf epscladdingsmf
f = k0.*asmf.*(sqrt((neff.^2) - epscladdingsmf));
end
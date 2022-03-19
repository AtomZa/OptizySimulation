function f = findwmmf(neff)
global k0 ammf epscladdingmmf
f = k0.*ammf.*(sqrt((neff.^2) - epscladdingmmf));
end
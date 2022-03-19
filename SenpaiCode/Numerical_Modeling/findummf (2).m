function f = findummf(neff)
global ammf k0 epscoremmf
f = k0.*ammf.*(sqrt(epscoremmf - (neff.^2)));
end

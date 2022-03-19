function f = characeqnsmf(neff)
u = findusmf(neff);
w = findwsmf(neff);
f = ((besselj(0,u))*(w)*(besselk(1,w))) - ((besselk(0,w))*(u)*(besselj(1,u)));  
end
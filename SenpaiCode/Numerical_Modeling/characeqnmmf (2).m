function f = characeqnmmf(neff)
u = findummf(neff);
w = findwmmf(neff);
f = ((besselj(0,u)).*(w).*(besselk(1,w,1))) - ((besselk(0,w,1)).*(u).*(besselj(1,u)));  
end

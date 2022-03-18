function f = fnnzno(vf)
global lambda nhost ncoremmf
f = ((1.781 + (0.059/(lambda^2)))*vf) + ((1 - vf)*nhost) - ncoremmf;
end
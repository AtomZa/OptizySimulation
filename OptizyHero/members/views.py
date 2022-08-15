from urllib import request
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from io import StringIO
from django.urls import reverse

import numpy as np
from matplotlib import pyplot as plt
from scipy.special import jv # Bessel function of the first kind
from scipy.special import kv # Modified Bessel function of the second kind
from scipy.special import kve # Exponentially scaled modified Bessel function of the second kind
from scipy.optimize import root_scalar # Root finding algorithm
import math

# Create your views here.
def index(request):
    # Data Analytics
    ncoresmf = 1.4504;
    ncladdingsmf = 1.447;
    asmf = 4.1;
    ncoremmf = 1.443;
    nmedium = 1.04641231840654;
    aair = 72.5;
    ammf = 62.5;
    zmmf = 57713;

    if request.method == "POST":
        ncoresmf = float(request.POST['ncoresmf']);
        ncladdingsmf = float(request.POST['ncladdingsmf']);
        asmf = float(request.POST['asmf']);
        ncoremmf = float(request.POST['ncoremmf']);
        nmedium = float(request.POST['nmedium']);
        aair = float(request.POST['aair']);
        ammf = float(request.POST['ammf']);
        zmmf = float(request.POST['zmmf']);

    ncladdingmmf = nmedium;
    epscoresmf = math.pow(ncoresmf, 2);
    epscladdingsmf = math.pow(ncladdingsmf, 2);
    epscoremmf = math.pow(ncoremmf, 2);
    epscladdingmmf = math.pow(ncladdingmmf, 2);
    lamb = 0.3;
    k0 = 2*math.pi/lamb;
    dx = 0.0001;
    ifinal = ((ncoresmf - ncladdingsmf)/dx);

    def characeqnsmf(neff,epscoresmf,epscladdingsmf,k0,asmf):
        return ((jv(0,k0*asmf*math.sqrt(epscoresmf - math.pow(neff,2))))*(k0*asmf*math.sqrt(math.pow(neff,2) - epscladdingsmf))*(kv(1,k0*asmf*math.sqrt(math.pow(neff,2) - epscladdingsmf)))) - ((kv(0,k0*asmf*math.sqrt(math.pow(neff,2) - epscladdingsmf)))*(k0*asmf*math.sqrt(epscoresmf - math.pow(neff,2)))*(jv(1,k0*asmf*math.sqrt(epscoresmf - math.pow(neff,2)))));

    neffarrsmf = 0;
    for i in range(round(ifinal)):
        try:
            sol = root_scalar(characeqnsmf, args=(epscoresmf, epscladdingsmf, k0, asmf), method='toms748', bracket=[ncoresmf-((i+1)*dx), ncoresmf-(i*dx)])
            neffarrsmf = sol.root
            break
        except:
            pass

    usmf =  k0*asmf*math.sqrt(math.pow(ncoresmf,2) - math.pow(neffarrsmf,2));   
    wsmf =  k0*asmf*math.sqrt(math.pow(neffarrsmf,2) - math.pow(ncladdingsmf,2));

    rarrsmf = np.arange(-asmf, asmf, 0.01);
    rarrcladding1 = np.arange(-ammf, -asmf-0.01, 0.01);
    rarrcladding2 = np.arange(asmf+0.01, ammf, 0.01);
    rarrair1 = np.arange(-aair, -ammf-0.01, 0.01);
    rarrair2 = np.arange(ammf+0.01, aair, 0.01);
    rarrall = np.concatenate((rarrair1, rarrcladding1, rarrsmf, rarrcladding2, rarrair2), axis=None)

    ecoresmf = jv(0,(usmf*rarrsmf/asmf));
    c = (jv(0,usmf))/(kv(0,wsmf));
    ecladdingsmf2 = c*(kv(0,(wsmf*rarrcladding2/asmf)));
    ecladdingsmf1 = np.zeros(np.size(ecladdingsmf2))
    eair1 = np.zeros(np.size(rarrair1));
    eair2 = np.zeros(np.size(rarrair2));
    eall =  np.concatenate((eair1, ecladdingsmf1, ecoresmf, ecladdingsmf2, eair2), axis=None);

    ifinal = ((ncoremmf - ncladdingmmf)/dx);
    neffarrmmf = np.zeros((1,round(ifinal)));

    def characeqnmmf(neff,epscoremmf,epscladdingmmf,k0,ammf):
        return ((jv(0,k0*ammf*math.sqrt(epscoremmf - math.pow(neff,2))))*(k0*ammf*math.sqrt(math.pow(neff,2) - epscladdingmmf))*(kve(1,k0*ammf*math.sqrt(math.pow(neff,2) - epscladdingmmf)))) - ((kve(0,k0*ammf*math.sqrt(math.pow(neff,2) - epscladdingmmf)))*(k0*ammf*math.sqrt(epscoremmf - math.pow(neff,2)))*(jv(1,k0*ammf*math.sqrt(epscoremmf - math.pow(neff,2)))));

    for i in range(round(ifinal)):
        try:
            sol = root_scalar(characeqnmmf, args=(epscoremmf, epscladdingmmf, k0, ammf), method='toms748', bracket=[ncladdingmmf + (i*dx), ncladdingmmf + ((i+1)*dx)])
            neffarrmmf[0,i] = sol.root
        except:
            neffarrmmf[0,i] = np.nan
            pass
    
    neffarrmmf= neffarrmmf[:, ~np.isnan(neffarrmmf).any(axis=0)]
    neffarrmmf = np.transpose(np.flip(neffarrmmf, axis=None))
    mmfmode = np.arange(0, np.size(neffarrmmf), 1)

    def return_graph():
        fig = plt.figure()
        plt.plot(rarrall,eall)

        imgdata = StringIO()
        fig.savefig(imgdata, format='svg')
        imgdata.seek(0)

        data = imgdata.getvalue()
        return data

    def mmf_graph():
        fig = plt.figure()
        plt.plot(mmfmode,neffarrmmf,'.')

        imgdata = StringIO()
        fig.savefig(imgdata, format='svg')
        imgdata.seek(0)

        data = imgdata.getvalue()
        return data

    # Web Application
    template = loader.get_template('myfirst.html')
    
    context = {
        'graph': return_graph(),
        'mmf': mmf_graph(),

        "ncoresmf": ncoresmf,
        "ncladdingsmf": ncladdingsmf,
        "asmf": asmf,
        "ncoremmf": ncoremmf,
        "nmedium": nmedium,
        "aair": aair,
        "ammf": ammf,
        "zmmf": zmmf,
    }
    return HttpResponse(template.render(context, request))
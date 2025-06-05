#!/usr/bin/env python


import sys
from scipy import special

# pice-wise power-law, e.g. Kroupa IMF
class c_imf_ppl:
    def __init__(self, edges, indeces):
        self.edges = edges
        self.indeces = indeces
        self.aconst = [0, ] * len(indeces)

        self.aconst[0] = 1.0
        for i in range(1,len(self.aconst)):
            self.aconst[i] = self.aconst[i-1]*self.edges[i]**(self.indeces[i-1]-self.indeces[i])
            #print '#', i, self.aconst[i]

        # fix alpha = -1 and -2 by a dirthy trick, 
        # because normalization coeffs diverge for these values
        for i in range(len(self.indeces)):
            if self.indeces[i] == -1.0 or self.indeces[i] == -2:
                self.indeces[i] += 1e-10

        # params for cumulative inverted IMF
        self.edgesCI = [0., ] * len(edges)
        #self.edgesCI[-1] = 1.0
        self.norm_number(1.0)
        for i in range(1, len(self.edges)):
            self.edgesCI[i] = self.edgesCI[i-1] + self.int_imf(self.edges[i-1], self.edges[i])



    def int_imf(self, m1, m2):
        intimf = 0
        for i in range(len(self.indeces)):
            if (m1 < self.edges[i+1]) and (m2 > self.edges[i]):
                mlo = max(m1, self.edges[i])
                mhi = min(m2, self.edges[i+1])
                intimf += self.aconst[i] * (1.0/(1.0+self.indeces[i])) \
                * (mhi**(1.0+self.indeces[i]) - mlo**(1.0+self.indeces[i]))
        return intimf

    def int_mass(self, m1, m2):
        intmass = 0
        for i in range(len(self.indeces)):
            if (m1 < self.edges[i+1]) and (m2 > self.edges[i]):
                mlo = max(m1, self.edges[i])
                mhi = min(m2, self.edges[i+1])
                intmass += self.aconst[i] * (1.0/(2.0+self.indeces[i])) \
                * (mhi**(2.0+self.indeces[i]) - mlo**(2.0+self.indeces[i]))
        return intmass

    def norm_number(self, n):
        norm = self.int_imf(self.edges[0], self.edges[-1])
        self.aconst = list(map(lambda x: x*n/norm, self.aconst))

    def norm_mass(self, mass):
        norm = self.int_mass(self.edges[0], self.edges[-1])
        self.aconst = list(map(lambda x: x*mass/norm, self.aconst))

    def cumInv(self, x):
        for i in range(len(self.indeces)):
            if (x < self.edgesCI[i+1]) and (x > self.edgesCI[i]):
                return (x*(self.edges[i+1]**(-self.indeces[i]+1.0)-self.edges[i]**(-self.indeces[i]+1.0)) \
                + self.edges[i]**(-self.indeces[i]+1.0))**(1.0/(-self.indeces[i]+1.0))

# Maschberger IMF - a continuous heavy-tailed log-normal distribution 
# with slopes of the tails similar to the Kroupa IMF
class c_imf_maschberger:
    def __init__(self, edges, indeces=(0.48, 2.3), mu=0.2):
        self.edges = edges
        self.alpha = indeces[1]
        self.gamma = indeces[0]
        self.beta = (self.gamma-self.alpha)/(1-self.alpha)
        self.mu = mu
        self.aconst = 1
        self.a = (2-self.alpha)/(1-self.alpha)
        self.b = self.beta - self.a

    def G(self, m):
        # an auxilliary function for computing the integral
        G = (1+(m/self.mu)**(1-self.alpha))**(1-self.beta)
        return G

    def t(self, m):
        # a transformation used for computing the mass integral
        t = ((m/self.mu)**(1-self.alpha))/(1+(m/self.mu)**(1-self.alpha))
        return t

    def int_imf(self, m1, m2):
        intimf = self.aconst*(self.mu/((1-self.alpha)*(1-self.beta)))*(self.G(m2)-self.G(m1))
        return intimf

    def int_mass(self, m1, m2):        
        intmass = self.aconst*(self.mu)**2/(1-self.alpha)*special.beta(self.a, self.b)*\
            (special.betainc(self.a, self.b, self.t(m2))-special.betainc(self.a, self.b, self.t(m1)))
        return intmass
    
    def norm_number(self, n):
        self.aconst = n/self.int_imf(self.edges[0], self.edges[1])

    def norm_mass(self, mass):
        self.aconst = mass/self.int_mass(self.edges[0], self.edges[1])

    def cumInv(self, x):
        if (x>0) and (x<1):
            return self.mu*((x*(self.G(self.edges[1])-self.G(self.edges[0]))+self.G(self.edges[0])) \
                **(1/(1-self.beta))-1)**(1/(1-self.alpha))

# parses cml arguments imfbin_str and alpha_str
# and sets list of imf mass bins 'imfbin' and power law coefficients 'alpha'
def parse_imfarg(imfbin_str, alpha_str, itype='ppl'):
    # convert imfbin_str to list imfbin
    comspl = imfbin_str.split(',')
    imfbin = list(map(lambda x: float(x), comspl))
    if itype.lower() == 'kroupa' or itype.lower() == 'salpeter' or itype.lower() == 'maschberger':
        if len(imfbin) != 2:
            print('Incorrect imfbin format.')
            print('Should be coma separated list.')
            print('Kroupa, Salpeter and Maschberger IMF can have only 2 bins: minimum and maximum mass.')
            sys.exit()

    if   itype.lower() == 'kroupa':
        imfbin = (imfbin[0], 0.08, 0.5, imfbin[1])
        alpha = (-0.3,-1.3,-2.3)
    elif itype.lower() == 'salpeter':
        alpha = (-2.35, )
    elif itype.lower() == 'maschberger':
        comspl = alpha_str.split(',')
        alpha = list(map(lambda x: float(x), comspl))
        if len(alpha) != 2:
            print('Only 2 alpha parameters can be set for the Maschberger IMF.')
    else:
        comspl = alpha_str.split(',')
        alpha = list(map(lambda x: float(x), comspl))
        if (len(alpha)+1) != len(imfbin):
            print('Mismatch between bin numbers of imfbin and alpha.')
            print('Number of imfbin elements must be by 1 greater than number of alpha elements')
            print(itype)
            sys.exit()

    return (imfbin, alpha)
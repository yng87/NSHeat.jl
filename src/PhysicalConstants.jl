#module PhysicalConstants
"""
Physical constatns and Unit conversion
"""
#export kB, hbarc, MeVToerg, MeVTofminv, mn, mp, me, mmu, G, c, Msun, sinWsq, gA, yrTosec, sigmaSB

const kB = 1.38e-16 #erg/K
const hbarc = 197.0 #MeV*fm
const MeVToerg = 1.6e-6 
const MeVTofminv = 1/hbarc #1/fm
const mn = 940. #MeV
const mp = 938. #MeV
const me = 0.511 #MeV
const mmu = 106. #MeV
const G = 6.671e-11 #m^3 kg^-1 s^-2
const c = 3.0e+8 #m/s
const Msun = 1.989e+30 #kg
const sinWsq= 0.231
const gA = 1.26
const yrTosec = 3.1536e+7 
const sigmaSB = 0.567 #erg s^-1 m^-2 K^-4

#end



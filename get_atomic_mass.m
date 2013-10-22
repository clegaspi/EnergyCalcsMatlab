function [ mass ] = get_atomic_mass( atomic_number )
% Returns atomic mass in atomic units

switch atomic_number
    case 1		% H, Hydrogen
        mass = 1.0079;
    case 2		% He, Helium
        mass = 4.0026;
    case 3		% Li, Lithium
        mass = 6.9410;
    case 4		% Be, Beryllium
        mass = 9.0122;
    case 5		% B, Boron
        mass = 10.8110;
    case 6		% C, Carbon
        mass = 12.0107;
    case 7		% N, Nitrogen
        mass = 14.0067;
    case 8		% O, Oxygen
        mass = 15.9994;
    case 9		% F, Fluorine
        mass = 18.9984;
    case 10		% Ne, Neon
        mass = 20.1797;
    case 11		% Na, Sodium
        mass = 22.9897;
    case 12		% Mg, Magnesium
        mass = 24.3050;
    case 13		% Al, Aluminum
        mass = 26.9815;
    case 14		% Si, Silicon
        mass = 28.0855;
    case 15		% P, Phosphorus
        mass = 30.9738;
    case 16		% S, Sulfur
        mass = 32.0650;
    case 17		% Cl, Chlorine
        mass = 35.4530;
    case 18		% Ar, Argon
        mass = 39.9480;
    case 19		% K, Potassium
        mass = 39.0983;
    case 20		% Ca, Calcium
        mass = 40.0780;
    case 21		% Sc, Scandium
        mass = 44.9559;
    case 22		% Ti, Titanium
        mass = 47.8670;
    case 23		% V, Vanadium
        mass = 50.9415;
    case 24		% Cr, Chromium
        mass = 51.9961;
    case 25		% Mn, Manganese
        mass = 54.9380;
    case 26		% Fe, Iron
        mass = 55.8450;
    case 27		% Co, Cobalt
        mass = 58.9332;
    case 28		% Ni, Nickel
        mass = 58.6934;
    case 29		% Cu, Copper
        mass = 63.5460;
    case 30		% Zn, Zinc
        mass = 65.3900;
    case 31		% Ga, Gallium
        mass = 69.7230;
    case 32		% Ge, Germanium
        mass = 72.6400;
    case 33		% As, Arsenic
        mass = 74.9216;
    case 34		% Se, Selenium
        mass = 78.9600;
    case 35		% Br, Bromine
        mass = 79.9040;
    case 36		% Kr, Krypton
        mass = 83.8000;
    case 37		% Rb, Rubidium
        mass = 85.4678;
    case 38		% Sr, Strontium
        mass = 87.6200;
    case 39		% Y, Yttrium
        mass = 88.9059;
    case 40		% Zr, Zirconium
        mass = 91.2240;
    case 41		% Nb, Niobium
        mass = 92.9064;
    case 42		% Mo, Molybdenum
        mass = 95.9400;
    case 43		% Tc, Technetium
        mass = 98.0000;
    case 44		% Ru, Ruthenium
        mass = 101.0700;
    case 45		% Rh, Rhodium
        mass = 102.9055;
    case 46		% Pd, Palladium
        mass = 106.4200;
    case 47		% Ag, Silver
        mass = 107.8682;
    case 48		% Cd, Cadmium
        mass = 112.4110;
    case 49		% In, Indium
        mass = 114.8180;
    case 50		% Sn, Tin
        mass = 118.7100;
    case 51		% Sb, Antimony
        mass = 121.7600;
    case 52		% Te, Tellurium
        mass = 127.6000;
    case 53		% I, Iodine
        mass = 126.9045;
    case 54		% Xe, Xenon
        mass = 131.2930;
    case 55		% Cs, Cesium
        mass = 132.9055;
    case 56		% Ba, Barium
        mass = 137.3270;
    case 57		% La, Lanthanum
        mass = 138.9055;
    case 58		% Ce, Cerium
        mass = 140.1160;
    case 59		% Pr, Praseodymium
        mass = 140.9077;
    case 60		% Nd, Neodymium
        mass = 144.2400;
    case 61		% Pm, Promethium
        mass = 145.0000;
    case 62		% Sm, Samarium
        mass = 150.3600;
    case 63		% Eu, Europium
        mass = 151.9640;
    case 64		% Gd, Gadolinium
        mass = 157.2500;
    case 65		% Tb, Terbium
        mass = 158.9253;
    case 66		% Dy, Dysprosium
        mass = 162.5000;
    case 67		% Ho, Holmium
        mass = 164.9303;
    case 68		% Er, Erbium
        mass = 167.2590;
    case 69		% Tm, Thulium
        mass = 168.9342;
    case 70		% Yb, Ytterbium
        mass = 173.0400;
    case 71		% Lu, Lutetium
        mass = 174.9670;
    case 72		% Hf, Hafnium
        mass = 178.4900;
    case 73		% Ta, Tantalum
        mass = 180.9479;
    case 74		% W, Tungsten
        mass = 183.8400;
    case 75		% Re, Rhenium
        mass = 186.2070;
    case 76		% Os, Osmium
        mass = 190.2300;
    case 77		% Ir, Iridium
        mass = 192.2170;
    case 78		% Pt, Platinum
        mass = 195.0780;
    case 79		% Au, Gold
        mass = 196.9665;
    case 80		% Hg, Mercury
        mass = 200.5900;
    case 81		% Tl, Thallium
        mass = 204.3833;
    case 82		% Pb, Lead
        mass = 207.2000;
    case 83		% Bi, Bismuth
        mass = 208.9804;
    case 84		% Po, Polonium
        mass = 209.0000;
    case 85		% At, Astatine
        mass = 210.0000;
    case 86		% Rn, Radon
        mass = 222.0000;
    case 87		% Fr, Francium
        mass = 223.0000;
    case 88		% Ra, Radium
        mass = 226.0000;
    case 89		% Ac, Actinium
        mass = 227.0000;
    case 90		% Th, Thorium
        mass = 232.0381;
    case 91		% Pa, Protactinium
        mass = 231.0359;
    case 92		% U, Uranium
        mass = 238.0289;
    case 93		% Np, Neptunium
        mass = 237.0000;
    case 94		% Pu, Plutonium
        mass = 244.0000;
    case 95		% Am, Americium
        mass = 243.0000;
    case 96		% Cm, Curium
        mass = 247.0000;
    case 97		% Bk, Berkelium
        mass = 247.0000;
    case 98		% Cf, Californium
        mass = 251.0000;
    case 99		% Es, Einsteinium
        mass = 252.0000;
    case 100		% Fm, Fermium
        mass = 257.0000;
    case 101		% Md, Mendelevium
        mass = 258.0000;
    case 102		% No, Nobelium
        mass = 259.0000;
    case 103		% Lr, Lawrencium
        mass = 262.0000;
    case 104		% Rf, Rutherfordium
        mass = 261.0000;
    case 105		% Db, Dubnium
        mass = 262.0000;
    case 106		% Sg, Seaborgium
        mass = 266.0000;
    case 107		% Bh, Bohrium
        mass = 264.0000;
    case 108		% Hs, Hassium
        mass = 277.0000;
    case 109		% Mt, Meitnerium
        mass = 268.0000;
    otherwise
        mass = 0.0;
end

end


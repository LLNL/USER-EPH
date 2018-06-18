#!/usr/bin/env octave

## unit conversion
cT = 10000; # temperature column is in 10^4K
cC = 6.2415e-7; # 10^5 J/m^3/K -> eV/A^3/K
cG = 6.2415e-7; # 10^17 W/m^3/K -> eV/ps/A^3/K

### read data from Zhigileis paper

## heat capacity
ce_in = load("Ce_Ni.dat");
N = length(ce_in(:,1));

## e-ph coupling
ge_in = load("G_Ni.dat");
M = length(ge_in(:,1));

#### OUTPUT RELATED
T = [0:100:1e+7]';
length(T);
T0 = ce_in(1,1)*cT;
T1 = ce_in(N,1)*cT;

ce0 = ce_in(1,2)*cC;
ce1 = ce_in(N,2)*cC;

ce = interp1(ce_in(:,1)*cT, ce_in(:,2)*cC, T, 'spline', 0.0);

ce = ce + ce0 * (T < T0);
ce = ce + ce1 * (T > T1);

out = [T, ce];
save('-ascii', 'C_e.data', 'out');

T0 = ge_in(1,1)*cT;
T1 = ge_in(M,1)*cT;

ge0 = ge_in(1,2)*cG;
ge1 = ge_in(M,2)*cG;

ge = interp1(ge_in(:,1)*cT, ge_in(:,2)*cG, T, 'spline', 0.0);

ge = ge + ge0 * (T < T0);
ge = ge + ge1 * (T > T1);

out = [T, ge];
save('-ascii', 'G_e.data', 'out');



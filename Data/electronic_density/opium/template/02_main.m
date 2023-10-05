#!/usr/bin/env octave

global interp_style = 'spline';

a0_to_A = 0.529177211;
ea03_to_eA3 = 6.748334491;

rcut = 5.5;
h = 0.9;

dr = 0.001;
r = [0.0:dr:rcut]';

XX_el = ZZ;

function fcut = fcut(r, rcut, h)
  rr = r < rcut; # if r > rc then fcut = 0.0
  rrr = ((r-rcut)/h).^4.0;
  rrrr = rrr + 1;
  rrrrr = rrr ./ rrrr;

  fcut = rrrrr .* rr;
endfunction

# this is in 4 pi * r^2 * rho(r)
core = load("xx_core_den.data");
valence = load("xx_valence_den.data");

core_rho = core(:,2) ./ core(:,1) ./ core(:,1) / 4.0 / pi * ea03_to_eA3;
core_r = core(:,1) * a0_to_A;

valence_rho = valence(:,2) ./ valence(:,1) ./  valence(:,1) / 4.0 / pi * ea03_to_eA3;
valence_r = valence(:,1) * a0_to_A;

val_rho = interp1(valence_r, valence_rho, r, interp_style, 0.0);
val_r = r;

all_rho = interp1(core_r, core_rho, r, interp_style, 0.0) + val_rho;
all_r = r;

val_rho = val_rho .* fcut(r, rcut, h);
all_rho = all_rho .* fcut(r, rcut, h);

# fix zeroth element
all_rho(1) = all_rho(2);

#plot(r, val_rho);
#pause();

#plot(r, all_rho);
#pause();

printf("Number of electrons in XX:      %.6f\n", XX_el)
N_el = 4.0 * pi * sum(r(:).^2 .* all_rho(:)) * dr;
printf("Integrated number of electrons: %.6f\n", N_el)

output = [r, all_rho * (XX_el/N_el)];
save("-ascii", "xx_rho.data", 'output');



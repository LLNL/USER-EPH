#!/usr/bin/env octave

### INPUT ###
kB = 8.617333262145e-5;

Zs = [29]; # use Zs = [1; 2; 3]; for multiple elements
Elems = ["Cu"]; # use Elems = ["A"; "B"; "C"] for multiple elements
N = length(Zs); # number of types
N_pairs = 1;

N_atoms = [2];

if(N > 1)
  N_pairs = (N+1)*(N-1) / 2;
endif

r_max = 5.0; # cutoff for spatial correlations
N_r = 1001; # number of points for table
dr = r_max / (N_r - 1);
r = [0.0:dr:r_max];

E_max = 1.0; # maximum energy in eV
N_E = 1001; # number of points for energy
dE = E_max / (N_E - 1);
E = [0.0:dE:E_max];
### END INPUT ###

### ARRAYS ###
rho_r = zeros(N, N_r); # spatial correlation
T_E = zeros(N, N_E); # electronic temperature
K_E = zeros(N_pairs, N_E);
### END ARRAYS ###

### INTERPOLATION ###
### RHO(R) ###
for i = [1:1:N]
  fn = sprintf("data/rho_%d.data", i);
  data = load(fn);
  r_in = data(:, 1);
  rho_in = data(:, 2);
  rho_interp = interp1(r_in, rho_in, r, "pchip", "extrap");
  rho_r(i, :) = rho_interp(:);
endfor
### END RHO_R ###

### T(E) ###
for i = [1:1:N]
  fn = sprintf("data/energy_%d.data", i);
  data = load(fn);
  T_in = data(:, 1);
  E_in = data(:, 2) / N_atoms(i);
  E_shift = E_in(1);

  T_E_interp = interp1(E_in - E_shift, T_in, E, "pchip", "extrap");
  T_E(i, :) = T_E_interp(:);
endfor
### END T(E) ###

### KAPPA(E) ###
k = 1;
for i = [1:1:N]
  for j = [i:1:N]
    fn = sprintf("data/kappa_%d_%d.data", i, j);
    data = load(fn);
    E_in = data(:, 1) / kB;
    K_in = data(:, 2);
    
    K_E_interp = interp1(E_in, K_in, E, "pchip", "extrap");
    K_E(k, :) = K_E_interp(:);

    k = k + 1;
  endfor
endfor
### END KAPPA(E) ###
### END INTERPOLATION ###

### SAVE ###
fout = fopen("out.kappa", "w");

fprintf(fout, "# comment line 1\n");
fprintf(fout, "# comment line 2\n");
fprintf(fout, "# comment line 3\n");

fprintf(fout, "%d ", N);

for i = [1:1:N]
  fprintf(fout, " %s", Elems(i, :));
endfor

fprintf(fout, "\n");
fprintf(fout, "%d %f %f %d %f %f\n", N_r, dr, r_max, N_E, dE, E_max);

for i = [1:1:N]
  fprintf(fout, "%d\n", Zs(i));
  
  for j = [1:1:N_r]
    fprintf(fout, "%f\n", rho_r(i, j));
  endfor
  fprintf(fout, "\n");

  for j = [1:1:N_E]
    fprintf(fout, "%f\n", T_E(i, j));
  endfor
  fprintf(fout, "\n");
endfor

for i = [1:1:N_pairs]
  for j = [1:1:N_E]
    fprintf(fout, "%f\n", K_E(i, j));
  endfor
  fprintf(fout, "\n");
endfor

fclose(fout);
### END SAVE ###


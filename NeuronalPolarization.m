% Matlab code for neural polarization
% Toriyama et al, Molecular Systems Biology, 2010
% DOI: 10.1038/msb.2010.51

% This script is the main code and uses the following sub scripts
% - script_model_setup.m
% - script_make_figure.m
%
%

clear variables;
close all;

%%%%%%%%%%%%%%%%%%%
% load setup
script_model_setup;

%%%%%%%%%%%%%%%%%%%
% main loop

for t = 2:t_end

  TIME = t .* dt;

  for k = 1:NEURITE_NUM

    % length dynamics
    Delta_A = fvL(L(t-1,k)) - fvC(C(t-1,k));
    DL = v_p .* M(t-1) - v_n .* exp(Delta_A);

    L(t,k) = L(t-1,k) + dt .* DL;

    % wave
    wave = W(t-1,k) .* C_soma(t-1);

    % concentration dynamics
    TRSPT = - ADV .* (C(t-1,k) - C_soma(t-1)) ./ L(t-1,k) + wave;
    C(t,k) = C(t-1,k) + dt .* TRSPT;

  end

  %%%%%
  Ec = hill(TIME,shootin_max,shootin_base,shootin_thalf,shootin_hill);
  Em = hill(TIME,tubulin_max,tubulin_base,tubulin_thalf,tubulin_hill);
  %%%%%

  % M
  M(t) = Em - sum(L(t,:)) * rho;

  % C
  C_soma(t) = Ec - sum(C(t,:))/soma_gc_volume_ratio;

end

%%%%%%%%%%%%%%%%%%%
% draw figures

script_make_figure;


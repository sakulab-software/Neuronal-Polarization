% Matlab code for neural polarization
% Toriyama et al, Molecular Systems Biology, 2010
% DOI: 10.1038/msb.2010.51

% This script is the sub code for
% - NeuronalPOlarization.m
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters

INPUT_NUM = 1000;
TIME_LIMIT = 8000;
NEURITE_NUM = 4;

dt = 0.1;
t_end = round(TIME_LIMIT/dt);
out_interval = 25;

%-----------------------------------------------------------
% Parameters in molecular expression formula

% conversion parameter between length and tubulin conc
rho = 1/300;  % relative concentration  (= 0.048/15)
% rho = 0.048;  % adjust to 15 uM

% Total shootin expression parameter
shootin_max   = 1.9;
shootin_base  = 0.07;
shootin_thalf = 3199;
shootin_hill  = 2.5;

% Total tubulin expression parameter
tubulin_max   = 2.9;       % relative concentration
% tubulin_max   = 2.9 *15;   % adjust to 15 uM
tubulin_base  = 1.1;       % relative concentration
% tubulin_base  = 1.1 *15;   % adjust to 15 uM
tubulin_thalf = 4480;
tubulin_hill  = 2.5;


%-----------------------------------------------------------
% Parameters in ODE

% dC/dt
ADV = 8.26;		                  % decay parameter of C at growth cone 
soma_gc_volume_ratio = 10;      % soma/growth cone
wave_gc_volume_ratio = 1;       % wave/growth cone

% dL/dt
v_p = 0.059 * 15;   % delta*kon  relative concentration
% v_p = 0.059;        % delta*kon  adjust to 15 uM
v_n = 0.68;    % delta*koff
% kappa * Fs
as  = 5.10;
Ks  = 2.15;
hc  = 4.30;
% kappa * Fl
al  = 2.39;
Kl  = 38.8;
L0  = 9.94;

%-----------------------------------------------------------
% Paramters for random wave

W_scale = 5.0;			% C(+) wave scale parameters
W_sigma = 0.4;			% wave width

% Wave intensity (gamma distrib)
intensity_a = 4.1921;
intensity_b = 0.2385;

% Inter wave interval (gamma distrib)
iwi_a = 2.2678;
iwi_b = 7.5054;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions

% Traction force function of Shootin1 conc
fvC = @(x)(as .* x .^ hc ./ (Ks .^ hc + x .^ hc));
% Tension function of neurite length
fvL = @(x)(al .* log(x ./L0) ./ (log(Kl ./ L0) + log(x ./ L0)));
% Sigmoid function
sig = @(t,a,b,c,baseline)( a ./ (1 + exp(-b .* (t-c))) + baseline);
% Hill formula
hill = @(t,a,b,k,h)(a .* t .^ h/(k .^ h + t .^h) + b);
% Wave shape (Gaussian)
wave_func = @(x)(exp(-(x ./ W_sigma) .^ 2 ./ 2.0) ./sqrt(2 .* pi) ./ W_sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables initialization

rng('shuffle');

% Neurite length
L      = zeros(t_end,NEURITE_NUM);
% Shootin1 conc in the growth cone
C      = zeros(t_end,NEURITE_NUM);
% Shootin1 conc at the cell body
C_soma = zeros(t_end,1);
% Transport by waves
W      = zeros(t_end,NEURITE_NUM);
% Free tubuline conc (uniform for cell body and growth cones)
M      = zeros(t_end,1);

% setup initial values
L(1,:)    = repmat(L0,1,NEURITE_NUM);
M(1)      = hill(dt,tubulin_max,tubulin_base,tubulin_thalf,tubulin_hill);

C_init    = hill(dt,shootin_max,shootin_base,shootin_thalf,shootin_hill);
C(1,:)    = repmat(C_init,1,NEURITE_NUM);
C_soma(1) = C_init;

% Calculation of transport time course (superposition of random waves)
A = tril(ones(INPUT_NUM,INPUT_NUM),0);
kernel = wave_func(-W_sigma*5:dt:W_sigma*5);
lag = round(W_sigma .*5 ./ dt);

for k = 1:NEURITE_NUM

  %----- TIMING
  iwis = gamrnd(iwi_a,iwi_b,INPUT_NUM,1);            % inter-wave intervals
  timings = A * iwis;                                % wave timings
  %-----

  timings_index = round(timings ./ dt);              % convert index
  timings_index(timings_index > t_end) = [];
  timings_index(timings_index <= W_sigma .*5 ./ dt) = [];
  a = zeros(t_end,1);
  a(timings_index) = 1;                              % pulse series(size 1)

  %----- INTENSITY
  % gamma
  amp = gamrnd(intensity_a, intensity_b, t_end, 1);  % intensity
  %-----

  amp = wave_gc_volume_ratio .* W_scale .* amp ./ (intensity_a .* intensity_b);
  a = a .* amp;                                      % weighted pulse
  W(:,k) = conv(a,kernel,'same');                    % gaussian convolitoin
end

clear A kernel amp iwis timings a;


% Matlab code for neural polarization
% Toriyama et al, Molecular Systems Biology, 2010
% DOI: 10.1038/msb.2010.51

% This script is the sub code for
% - NeuronalPOlarization.m
%
%

before = 1000;   % Interval to display BEFORE polarization time (minutes)
after  = 1500;   % Interval to display AFTER polarization time (minutes)

%-------------------------------------------------
% setup 

f = figure(1); clf;
f.Position = [1050 150 750 820]; 
set(gcf, 'Color', 'w');

corder = [[1 0 0];
          [0 0 1];
          [0 0.5 0];
          [1 0.75 0];
          [0.5 0.3 0.2];
          [0.8 0 0.8];
          [0.25 0.25 0.25]
         ];  
set(0, 'defaultAxesColorOrder', corder);

%-------------------------------------------------

% Find the longest
L_last = L(end,:);
[~, idx] = sort(L_last, 'descend');
L = L(:,idx);
C = C(:,idx);

[i, ~] = ind2sub(size(L),find(L >200 & L<202, 1 ));
dt = 0.1;
from = i*dt - before; toto = min(from + after,TIME_LIMIT);  

t = dt:dt:TIME_LIMIT;
r = [round(from/dt):round(toto/dt)];

legend_entries = cell(1, NEURITE_NUM);

%---
subplot(3,1,1);
hold on;
for i = 1:NEURITE_NUM
  plot(t(r), L(r,i));
  legend_entries{i} = ['Neurite ' num2str(i)];
end
legend(legend_entries, 'Location','northwest');

axis([from toto 0 200]);
xlabel('Time (min)');
ylabel('Length (um)');
title('Neurite lengths');

%---
subplot(3,1,2);
hold on;
for i = 1:NEURITE_NUM
  plot(t(r), C(r,i));
  legend_entries{i} = ['Neurite ' num2str(i)];
end
legend(legend_entries, 'Location','northwest');

axis([from toto 0 15]);
xlabel('Time (min)');
ylabel('Shootin1 concentration');
title('Shootin1 relative concentrations in growth cones');

%---
subplot(3,1,3);
hold on;
yyaxis left;
plot(t(r), C_soma(r));
axis([from toto 0 3]);
ylabel('Shootin1 concentration');
ax = gca;
ax.YColor = 'r';

yyaxis right;
plot(t(r), M(r), 'k-');
axis([from toto 0 max(M)*1.5]);
ylabel('Tubuline concentration');
ax = gca;
ax.YColor = 'k';

xlabel('Time (min)');
legend('Shootin1 in cell body', 'Free tubulin', 'Location','northwest');
title('Shootin1 and tubulin');




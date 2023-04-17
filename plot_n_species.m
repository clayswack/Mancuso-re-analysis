%% plot_n_species

%% Objective
% Load an output from the mancuso model evaluated in jupyter lab
% The output is a 4D array (dilution, freq, simulation_number, resource or bacteria population)

%% Actually, forget about it. 
% For now, calculate shannon in python. 
% Import the values themselves into matlab

%% Setup
clear
close all

%% Get data
% load test_20230410_1.mat % the full 4D array

d = 0:0.1:1.1; % dilution
sha = getSha("C:\Users\clayswack\OneDrive\OneDrive Docs\Purdue\2023 Bioinformatics class\Project analysis\Modeling\nASVs_0411_1.txt", [1, Inf]);

% Sha is 12 by 7
% 12 dilution rates by 7 frequencies

%% Make plot
markers = {'o', 's', 'd', '*', '+', 'pentagram', '>'}
figure

for i = 1:7
    plot(d, sha{:,i}, markers{i}, 'LineWidth', 1.5, 'LineStyle', '-')
    hold on
end

%% Plot style
ylabel('Number species with N_i > 0.05 ', 'FontSize', 14)
xlabel('Mean dilution rate (hr^{-1})', 'FontSize', 14)
 
ax = gca;
ax.YRuler.TickLabelFormat = '%.1f'; %enforce consistent xtick and ytick decimals
ax.XRuler.TickLabelFormat = '%.1f';

legTxt = {'1 dilution / day', '2 dilution / day', '4 dilution / day',...
    '8 dilution / day', '16 dilution / day', '32 dilution / day', '64 dilution / day'}
leg = legend(legTxt);
leg.Location = "eastoutside";
leg.FontSize = 12;

set(gcf, 'Position', [744	630	673	420])

set(gca,'linewidth',1.25) %change box line width











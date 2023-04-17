%% plot_alpha_rafefaction

%% Objective:
% Graph shannon diversity vs alpha rarefaction depth
% For this graph, break it out by dilution rate and frequency separately

%% Setup
clear
close all
clc

%% Get data
data = importfile2("C:\Users\clayswack\OneDrive\OneDrive Docs\Purdue\2023 Bioinformatics class\Figs\Rarefaction fig\alpha_rarefaction_shannon_matlab.xlsx", "alpha_rarefaction_shannon_matla", [2, 630]);
head(data)

%% Average all iterations
for i = 1:height(data)
    data.iterAvs(i) = mean(data{i,5:end});
end

dataSparse = table;
dataSparse.dilution = data.dilution;
dataSparse.freq = data.freq;
dataSparse.rep = data.rep;
dataSparse.depth = data.depth;
dataSparse.y = data.iterAvs;

head(dataSparse)

%% Average the 4 replicates at each condition
avsReps = grpstats(dataSparse, {'dilution', 'freq', 'depth'}, {'mean'}, 'DataVars', {'y'},...
    'VarNames', {'dilution', 'freq', 'depth', 'grpCount', 'y'})

% errsReps = grpstats(dataSparse, {'dilution', 'freq', 'depth'}, {'std'}, 'DataVars', {'y'},...
%     'VarNames', {'dilution', 'freq', 'depth', 'grpCount', 'err'});
% avsReps.err = errsReps.err;

%% Next average out the different frequencies to plot by dilution
avsRepsFreqs = grpstats(avsReps, {'dilution', 'depth'}, {'mean'}, 'DataVars', {'y'},...
    'VarNames', {'dilution', 'depth', 'grpCount', 'y'})

%% Make plot by dilution
d1 = avsRepsFreqs(avsRepsFreqs.dilution == 0.1,:);
d2 = avsRepsFreqs(avsRepsFreqs.dilution == 0.2,:);
d3 = avsRepsFreqs(avsRepsFreqs.dilution == 0.3,:);
d4 = avsRepsFreqs(avsRepsFreqs.dilution == 0.4,:);

subplot(1,2,1)
plot(d1.depth, d1.y, ':o', 'LineWidth', 1.5)
hold on
plot(d2.depth, d2.y, ':s', 'LineWidth', 1.5)
plot(d4.depth, d3.y, ':p', 'LineWidth', 1.5)
plot(d1.depth, d1.y, ':pentagram', 'LineWidth', 1.5)

% % Plot the value that was chosen
% plot([6840, 6840], [0, 1.2], 'k:', 'LineWidth', 3)

leg = legend('Dilution = 0.1/hr', 'Dilution = 0.2/hr', 'Dilution = 0.3/hr', 'Dilution = 0.4/hr');
leg.Location = 'SouthEast';
leg.FontSize = 12;

xlabel('Rarefaction depth (number of bases)', 'FontSize', 12)
ylabel('Shannon diversity', 'FontSize', 12)

ax = gca;
ax.YRuler.TickLabelFormat = '%.1f'; %enforce consistent xtick and ytick decimals

set(gca,'linewidth',1.25) %change box line width

%% Next average out the different dilutions to plot by frequency
avsRepsDilutions = grpstats(avsReps, {'freq', 'depth'}, {'mean'}, 'DataVars', {'y'},...
    'VarNames', {'freq', 'depth', 'grpCount', 'y'})

%% Now break it down by frequency
f1 = avsRepsDilutions(avsRepsDilutions.freq == 1,:);
f2 = avsRepsDilutions(avsRepsDilutions.freq == 4,:);
f3 = avsRepsDilutions(avsRepsDilutions.freq == 16,:);
f4 = avsRepsDilutions(avsRepsDilutions.freq == 600,:);

subplot(1,2,2)
plot(f1.depth, f1.y, ':o', 'LineWidth', 1.5)
hold on
plot(f2.depth, f2.y, ':s', 'LineWidth', 1.5)
plot(f3.depth, f3.y, ':p', 'LineWidth', 1.5)
plot(f4.depth, f4.y, ':pentagram', 'LineWidth', 1.5)

% % Plot the value that was chosen
% plot([6840, 6840], [0, 2.5], 'k:', 'LineWidth', 3)

leg = legend('Frequency = 1/day', 'Frequency = 4/day', 'Frequency = 16/day', 'Frequency = 600/day');
leg.Location = 'Best';
leg.FontSize = 12;

xlabel('Rarefaction depth (number of bases)', 'FontSize', 12)
ylabel('Shannon diversity', 'FontSize', 12)

ax = gca;
ax.YRuler.TickLabelFormat = '%.1f'; %enforce consistent xtick and ytick decimals

set(gca,'linewidth',1.25) %change box line width

%}

%% Overall plot style
set(gcf, 'Position', [744	630	994.600000000000	377.400000000000])













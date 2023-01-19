clear all, close all, clc

load('source_signal.mat');

%% settings
linewidth=1;
fs_label=12;
fs_tick=10;
fs_legend=10;
figsize=[100,100,400,250];

figure('color','w','position',figsize);
plot(h_star,'k');
set(gca,'FontSize',fs_tick);
xlim([0 512])
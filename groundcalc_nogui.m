clear all;
close all;
clc;

addpath('functions');

fname='Drawing1.txt';
delta=.5;
buffer=10;
cond_radius=0.0067056;
rho_top=1000;
rho_bottom=100;
h_top=1;
i_src=complex(3000,0);
t_fault=1000;
rho_cover=3000;
h_cover=0 / 100;
plot_surf=true;
plot_curr=true;
plot_touch=true;
plot_step=true;
plot_grid=true;

try
    runGroundCalcSolver(fname, delta, buffer, cond_radius, rho_top, rho_bottom, h_top, i_src, t_fault, rho_cover, h_cover, plot_surf, plot_curr, plot_touch, plot_step, plot_grid);
catch err
end
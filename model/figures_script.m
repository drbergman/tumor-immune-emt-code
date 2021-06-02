close all; clear all;
% This should be the script that built most of the figures
[~,~,~,keep_warmup] = standard_start();

new_pars.E_NK_up = 1.2;
new_pars.E_CTL_up = 3;
new_pars.E_Treg_up = 10;

new_pars.K3 = 200;
new_pars.TGFB_max = 700;
new_pars.E_NK_low = 1;
new_pars.E_CTL_low = 100;

new_pars.INFL_LOW_duration = 30;
new_pars.INFL_HIGH_duration = 30;

%% Figure 1B
N = 1;
T = 2000;
tracker_index = [1,2,3,4,6,7,10];

trial{1} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{1}.name = 'fig1b';
%% Figure 1C
tracker_index = [];
T = 2000;
N = 500;

trial{2} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{2}.name = 'fig1c';
%% Figure 3A, D

N = 1000;
T = 2000;
tracker_index = [1];

new_pars.mie = linspace(.2,.8,4);

trial{3} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{3}.name = 'fig3a';
%% Figure 3B, E
[tracker_index,T,N,keep_warmup] = standard_start();

N = 1000;
T = 2000;
tracker_index = [];

new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;

new_pars.INFL_LOW_duration = 30;
new_pars.INFL_HIGH_duration = 60;

new_pars.mie = .7;
new_pars.mga = linspace(.1,.4,4);

new_pars.E_NK_up = 1.2;
new_pars.E_CTL_up = 3;
new_pars.E_Treg_up = 10;

trial{4} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{4}.name = 'fig3b';
%% Figure 3C, F
[tracker_index,T,N,keep_warmup] = standard_start();

N = 1000;
T = 2000;
tracker_index = [];

new_pars.E_NK_low = 10;
new_pars.K2 = 2.5;
new_pars.E_CTL_low = 200;
new_pars.mie = .5;
new_pars.mga = .2;

new_pars.NKrec_up = 3;
new_pars.CTLrec_up = 3;
new_pars.Tregrec_up = 5;

new_pars.TGFB_MUT = .5;
new_pars.TGFB_Treg = linspace(.1,1,4);

trial{5} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{5}.name = 'fig3c';
%% Figure 4A, B
[tracker_index,T,N,keep_warmup] = standard_start();

N = 200;
T = 2000;
tracker_index = [];

new_pars.mie = linspace(.2,.8,4);
new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;
new_pars.mga = 0;
new_pars.INFL_LOW_duration = [0,7];
new_pars.INFL_HIGH_duration = [0,7];

trial{6} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{6}.name = 'fig4a';
%% Figure 4C, D
[tracker_index,T,N,keep_warmup] = standard_start();

N = 200;
T = 2000;
tracker_index = [];

new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;
new_pars.INFL_LOW_duration = [0,7];
new_pars.INFL_HIGH_duration = [0,7];

new_pars.mie = .9;
new_pars.mga = linspace(.1,.7,4);

trial{7} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{7}.name = 'fig4c';
%% Figure 5
[tracker_index,T,N,keep_warmup] = standard_start();

N = 200;
T = 2000;
tracker_index = [];

new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;

new_pars.INFL_LOW_duration = 30;
new_pars.INFL_HIGH_duration = 60;

new_pars.mie = linspace(.4,.8,11);
new_pars.mga = linspace(.1,.4,11);

new_pars.E_NK_up = 1.2;
new_pars.E_CTL_up = 3;
new_pars.E_Treg_up = 10;

new_pars.K3 = 200;
new_pars.TGFB_max = 700;

trial{8} = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
trial{8}.name = 'fig5';

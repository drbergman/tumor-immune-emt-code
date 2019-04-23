%% Figure 1B
N = 1;
T = 2000;
tracker_index = [1,2,3,4,6,7,10];

new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;

new_pars.mie = .9;
new_pars.mga = .2;

new_pars.K3 = 200;
new_pars.TGFB_max = 700;

trial.1b = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 1C
N = 500;

new_pars.mie = linspace(.4,.9,4);
new_pars.mga = linspace(.1,.4,4);

trial.1c = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 3A, D
N = 1000;
T = 2000;
tracker_index = [1,6];

new_pars.mie = linspace(.2,.8,4);
new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;
new_pars.mga = 0;

trial.3a = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 3B, E
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

trial.3b = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 3C, F
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

trial.3c = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 4A, B
N = 200;
T = 2000;
tracker_index = [];

new_pars.mie = linspace(.2,.8,4);
new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;
new_pars.mga = 0;
new_pars.INFL_LOW_duration = [0,7];
new_pars.INFL_HIGH_duration = [0,7];

trial.4a = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 4C, D
N = 200;
T = 2000;
tracker_index = [];

new_pars.E_NK_low = 10;
new_pars.E_CTL_low = 200;
new_pars.INFL_LOW_duration = [0,7];
new_pars.INFL_HIGH_duration = [0,7];

new_pars.mie = .9;
new_pars.mga = linspace(.1,.7,4);

trial.4c = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);

%% Figure 5
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

trial.5 = cohort_fn(new_pars,T,N,tracker_index,keep_warmup);
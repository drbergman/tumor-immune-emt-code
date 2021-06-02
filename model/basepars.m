function pars = basepars()

pars.p = .28;
pars.d_C = .14;
pars.mie = .6;
pars.mga = .2;

pars.apop_down = .3;
pars.immune_evasion = .48;
pars.prolif_up = .36;

pars.K0 = 80;
pars.K1 = 8;
pars.K2_low = 5;
pars.K3 = 200;
pars.K4 = 50;

pars.E_NK_low = .2;
pars.E_CTL_low = 4;
pars.E_NK_nonmut = 0;
pars.E_CTL_nonmut = .1;

pars.E_NK_up = 3;
pars.E_CTL_up = 3;
pars.E_Treg_up = 200;

pars.sigma_NK_low = 1.3;
pars.sigma_CTL_low = 100;
pars.sigma_Treg_low = 200;

pars.sigma_NK_up = 1;
pars.sigma_CTL_up = 1;
pars.sigma_Treg_up = 1;

pars.d_NK = .13;
pars.d_CTL = pars.d_NK/5;
pars.d_Treg = pars.d_CTL;

pars.k_EMT = .01;
pars.sigma = 6;

pars.TGFB_max = 500;
pars.TGFB_MUT = 5e-2;
pars.TGFB_Treg = 5e-1;

pars.RP_Cancer_line = .5;

pars.INFL_HIGH_duration = 30;
pars.INFL_LOW_duration = 30;

pars.Mes_threshold = .7;

pars.p_mutation_start = 1e-2;
pars.p_mutation_on = 1e-4;

pars.Warmup = 1000;
pars.Nmax = 400;
pars.t_step = 1;

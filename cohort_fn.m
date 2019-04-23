function trial = cohort_fn(new_pars,T,N,TI,KWU)

trial.N = N;
trial.T = T;
trial.tracker_index = TI;
trial.keep_warmup = KWU;

[tracked_names,trial.per_Nums] = track_ind2names(TI);

[trial.pars,trial.sub_cohort_sz] = buildpars(new_pars);

trial.pat = sub_cohort(T,trial.pars,...
            KWU,trial.sub_cohort_sz,N,tracked_names);
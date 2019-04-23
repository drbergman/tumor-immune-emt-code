function EMT = update_emt(EMT,TGFB,N,pars)

TGFB_received = pars.TGFB_max*(TGFB/pars.K3/(1+TGFB/pars.K3))/N+pars.sigma*randn(N,1)+(EMT-2);
EMTed = TGFB_received > 0;
METed = TGFB_received < 0;
EMT(EMTed) = 1-(1-EMT(EMTed)).*exp(-pars.k_EMT*(TGFB_received(EMTed)));
EMT(METed) = EMT(METed).*exp(2*pars.k_EMT*(TGFB_received(METed)));

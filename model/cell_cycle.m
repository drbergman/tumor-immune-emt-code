function [C,alive,Num,Mut,NK,CTL,Treg,TGFB,EMT,p_cell_mutate] = ...
    cell_cycle(C,alive,pars,Num,NK,CTL,Treg,TGFB,EMT,p_cell_mutate)

%% Build Probability Matrix
    p_apoptosis = pars.d_C*(1-pars.apop_down*(C(:,1)~=0));
    p_NK_clearance = (any(C,2).*((NK/(Num/pars.K1+NK))*pars.E_NK*(1/(1+Treg/pars.K2)))).*...
        (1-pars.mie*(EMT>pars.Mes_threshold)).*(1-pars.immune_evasion*(C(:,3)~=0));
    p_CTL_clearance = (any(C,2).*((CTL/(Num/pars.K1+CTL))*pars.E_CTL*(1/(1+Treg/pars.K2)))).*...
        (1-pars.mie*(EMT>pars.Mes_threshold)).*(1-pars.immune_evasion*(C(:,3)~=0));
    p_proliferation = pars.p*(1+pars.prolif_up*(C(:,2)~=0)).*...
        (1-pars.mga*(EMT>pars.Mes_threshold))*pars.K0/(pars.K0+Num);
    p_quiescence = ones(pars.Nmax,1)+...
        ...% This next line takes the lost prolif prob of Mes cells and adds
        ...% it to their prob to be quiescent
        pars.p*(1+pars.prolif_up*(C(:,2)~=0)).*(pars.mga*(EMT>pars.Mes_threshold))*pars.K0/(pars.K0+Num);
    
    p_matrix = [p_apoptosis,p_NK_clearance,p_CTL_clearance,p_proliferation,...
        p_quiescence];
    p_matrix = p_matrix./sum(p_matrix,2);
    
%% Cell apoptosis
    r = rand(pars.Nmax,1) - p_matrix(:,1);
    alive((r<0) & (alive)) = false;
    
%% Immune clearance
    r = r - p_matrix(:,2);
    NK_lysis = find((r<0) & (alive));
    alive(NK_lysis) = false;
    r = r - p_matrix(:,3);
    CTL_lysis = find((r<0) & (alive));
    alive(CTL_lysis) = false;
    
%% Cell proliferation
    r = r - p_matrix(:,4);
    proliferate = find((r<0) & (alive));
    
    NK = max(0,NK-length(NK_lysis));
    CTL = max(0,CTL-length(CTL_lysis));

%% Mutation
    if ~isempty(proliferate)
        num_proliferate = length(proliferate);
        mutate = rand(num_proliferate,1)<p_cell_mutate(proliferate);
        if ~isempty(find(mutate,1))
            C(sub2ind([pars.Nmax,3],proliferate(mutate),randi(3,size(find(mutate))))) = 1;
            p_cell_mutate(proliferate(mutate)) = 0;
        end
        p_cell_mutate(proliferate(~mutate)) = p_cell_mutate(proliferate(~mutate))+pars.p_mutation;
        new_cell_ind = find(~alive,num_proliferate);
        alive(new_cell_ind) = true;
        try
            C(new_cell_ind,:) = C(proliferate,:);
            EMT(new_cell_ind) = EMT(proliferate);
        catch
            error('Not enough space to put new cells');
        end
    end
    
    Num = sum(alive);
    Mut = sum(any(C,2) & alive);
    C(~alive,:) = 0;
    
%% EMT update
    EMT(alive) = update_emt(EMT(alive),TGFB,Num,pars);
    TGFB = Mut*pars.TGFB_MUT+Treg*pars.TGFB_Treg;
    
%% Immune system update
    NK = (NK-pars.sigma_NK/pars.d_NK)*exp(-pars.d_NK*pars.t_step)+pars.sigma_NK/pars.d_NK;
    activated_DC = length(NK_lysis)+length(CTL_lysis);
    CTL = (CTL-pars.sigma_CTL*activated_DC)*exp(-pars.d_CTL*pars.t_step)+...
        pars.sigma_CTL*activated_DC;
    treg_recruit = pars.sigma_Treg*TGFB/(1+TGFB/pars.K4);
    Treg = (Treg-treg_recruit*activated_DC)*exp(-pars.d_Treg*pars.t_step)+...
        treg_recruit*activated_DC;
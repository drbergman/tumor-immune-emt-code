function pat = patient(T,pars,KWU,TN)

pat.T2C = [];
pat.num_at_cancer = NaN;

cycle_length = max(1,pars.INFL_HIGH_duration+pars.INFL_LOW_duration);
cycle_shift = randi(cycle_length);

%% Mutational parameters
pars.p_mutation = 0;

%% Start with N cells, each healthy (0) or mutated (1).
Num = 100;

C = zeros(pars.Nmax,3);

alive = [true(Num,1);false(pars.Nmax-Num,1)];
Mut = sum(alive & any(C,2));
if Mut > 1
    disp(C(1,:))
end

%% Initialize an immune system with no adaptive immune cells
NK = 10;
CTL = 0;
Treg = 0;
TGFB = 0;

%% EMT parameters
EMT = [.1*rand(Num,1);zeros(pars.Nmax-Num,1)];

%% Mutational parameters
p_cell_mutate = zeros(pars.Nmax,1);

Warmup = pars.Warmup;

Mut_types = sum(C(alive,:)*[1;2;4]==(0:7))';
EMT_val_quarts = prctile(EMT(alive),[0;25;50;75;100]);
EMT_bincounts = histcounts(EMT(alive),(0:.1:1))';
Mes = sum(alive & (EMT>pars.Mes_threshold));

%% Values to track
for ti = 1:length(TN)
    pat.(TN{ti}) = zeros(numel(eval(TN{ti})),KWU*(Warmup+1)+(T+1));
    pat.(TN{ti})(:,1) = eval(TN{ti});
end

pars.E_NK = pars.E_NK_low;
pars.K2 = pars.K2_low;
pars.E_CTL = pars.E_CTL_low;
pars.sigma_NK = pars.sigma_NK_low;
pars.sigma_CTL = pars.sigma_CTL_low;
pars.sigma_Treg = pars.sigma_Treg_low;
%% Run T Cell Cycles
for i = -Warmup:T
    if (mod(i+cycle_shift,cycle_length) == 0) && (pars.INFL_HIGH_duration ~= 0) % Now INFL is high
        pars.E_NK = pars.E_NK_low*pars.NKeff_up;
        pars.E_CTL = pars.E_CTL_low*pars.CTLeff_up;
        pars.K2 = pars.K2_low/pars.Tregeff_up;
        pars.sigma_NK = pars.sigma_NK_low*pars.sigma_NK_up;
        pars.sigma_CTL = pars.sigma_CTL_low*pars.sigma_CTL_up;
        pars.sigma_Treg = pars.sigma_Treg_low*pars.sigma_Treg_up;
    elseif mod(i+cycle_shift,cycle_length) == pars.INFL_HIGH_duration % Now INFL is low
        pars.E_NK = pars.E_NK_low;
        pars.E_CTL = pars.E_CTL_low;
        pars.K2 = pars.K2_low;
        pars.sigma_NK = pars.sigma_NK_low;
        pars.sigma_CTL = pars.sigma_CTL_low;
        pars.sigma_Treg = pars.sigma_Treg_low;
    end

    if i == 1
        pars.p_mutation = pars.p_mutation_on;
        p_cell_mutate(alive) = (pars.p_mutation_start)*rand(Num,1);
    end
    
    [C,alive,Num,Mut,NK,CTL,Treg,TGFB,EMT,p_cell_mutate] = ...
        cell_cycle(C,alive,pars,Num,NK,CTL,Treg,TGFB,EMT,p_cell_mutate);
    
%% Update tracked values
    Mut_types = sum(C(alive,:)*[1;2;4]==(0:7))';
    EMT_val_quarts = prctile(EMT(alive),[0;25;50;75;100]);
    EMT_bincounts = histcounts(EMT(alive),(0:.1:1))';
    Mes = sum(alive & (EMT>pars.Mes_threshold));
    
    for ti = 1:length(TN)
        pat.(TN{ti})(:,KWU*(Warmup+1)+i+1) = eval(TN{ti});
    end
    
%% Check if cancer can be detected
    if (((Mut/Num >= pars.RP_Cancer_line) && Num >= 0) || Num == 0) && isempty(pat.T2C)
        
        pat.T2C = i;
        pat.num_at_cancer = Num;
        if isempty(TN)
            break
        end
    end
end

if isempty(pat.T2C)
    pat.T2C = Inf;
end

% rearrange so single mutations are in rows 2-4 and double mutations in 5-7
if isfield(pat,'Mut_types')
    pat.Mut_types([4,5],:) = pat.Mut_types([5,4],:);
end

% rearrange so quartiles go from top to bottom for graphing purposes
if isfield(pat,'EMT_val_quarts')
    pat.EMT_val_quarts = pat.EMT_val_quarts(end:-1:1,:);
end

function pat = sub_cohort(T,pars,keep_warmup,sub_cohort_sz,N,TN)
        
if ~isempty(sub_cohort_sz)
    colons = repmat({':'},1,length(sub_cohort_sz));
    P = pars;
    ind = find(structfun(@numel,pars)>1,1);
    names = fieldnames(pars);
    for i = sub_cohort_sz(1):-1:1
        P.(names{ind}) = pars.(names{ind})(i);
        pat(i,colons{:}) = sub_cohort(T,P,...
            keep_warmup,sub_cohort_sz(2:end),N,TN);
    end
else
    for j = 1:N
        pat(j) = patient(T,pars,keep_warmup,TN);
    end
end
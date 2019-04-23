function [pars,scz] = buildpars(new_pars)

pars = basepars();

names = fieldnames(new_pars);
for i = 1:length(names)
    if ~isfield(pars,names{i})
        warning('The new_pars field %s is not a field in basepars().',names{i})
    end
    pars.(names{i}) = new_pars.(names{i});
end

scz = structfun(@numel,pars);
scz = scz(scz>1);

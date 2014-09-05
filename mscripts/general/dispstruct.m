

function dispstruct(fpt,opt)


F = fieldnames(opt);
for i = 1:length(F)
  f = F{i};
  v = getfield(opt,{1},F{i});
  v = num2str(v);
  fprintf(fpt,'%24s:  %s\n',f,v);
end




function x = rinvgam(Shape,Scale)
%% Parametrize in terms of shape and scale %%
%% rather than shape and rate              %%


%% randg is in the statistics toolbox
if (nargin~=2)
  error('rinvgam:NotEnoughInputs','Wrong number of parameters.')
end
if ~(isscalar(Shape) && isscalar(Scale))
  error('rinvgam:WrongUsage','usage:rinvgam(Shape,Scale).')
end

%%x = Scale./randg(Shape,n,1);
x = Scale./randraw('gamma',[Shape],1);

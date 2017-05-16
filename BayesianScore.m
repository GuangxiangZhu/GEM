function [ score ] = BayesianScore( structure,proportions,C1,v0 )

% v0: The real volume
% v: The estimated volume of the computed structures
v=Volume(structure,proportions);
score=(1-C1)*v0/abs(v-v0);

end


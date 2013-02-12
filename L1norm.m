function D = L1norm(prob,prob2)

%% compute distance based on L1 norm divided by two to approximate the earth movers
% prob,prob2 : Two distributions

if min(min(prob),min(prob2)) < 0
    warning('Either or both distributions contain negative values');
end

if (length(prob)~= length(prob2))
    warning('Distributions do not have the same support');
end

D = sum(abs(prob - prob2))/2;
if D < 1e-15
    D = 0;
end
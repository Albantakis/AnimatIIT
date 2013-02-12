function ind = convi(M_pb)

%% M_pb is based on binary counting
N = length(M_pb);

ind = 1;
for i = 1:N
    j = M_pb(i);
    ind = ind + 2^(j-1);
end

end
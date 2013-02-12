function index = subsystem2index(subsystem)

N = length(subsystem);

index = 0;
for i = 1:N
    j = subsystem(i);
    index = index + 2^(j-1);
end

end
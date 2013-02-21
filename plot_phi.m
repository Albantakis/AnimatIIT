function [] = plot_phi(phi_vec,MIP_cell, M)

N = length(M);

string = cell(2^N-1,1);
phi_vec2 = zeros(2^N-1,1);
k = 1;
l = 1;
for i=1: N
    C = nchoosek(M,i);
    N_C = size(C,1);
    for j=1: N_C
        x = C(j,:);
        string{k} = num2str(x);
        if l <= size(MIP_cell,1)
            MIP = MIP_cell{l};
            y = combine(MIP{1,2},MIP{2,2});
            if  strcmp(num2str(x),num2str(y)) == 1
                phi_vec2(k) = phi_vec(l);
                l = l +1;
            end
        end
        k = k + 1;
    end
end

h = bar(1:2^N-1,phi_vec2);
set(h,'facecolor','black');
set(gca,'XTick',1:2^N-1)
set(gca,'XTickLabel',string) % uncomment this to have a binary valued x-axis
% rotateXLabels( gca(), 90)
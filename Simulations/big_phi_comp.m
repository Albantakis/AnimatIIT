function [Big_phi phi prob_cell] = big_phi_comp(M,x1,p,b_table,options)

%%  compute big phi in a subset M
% M: a subset of the whole system (can be the whole system itself)
% x1: given data about the current state
% x0: past state x1: current state
% p: transition probability matrix in the whole system (p(x1|x0))

N = length(M);
prob_cell = cell(3,1);

op_fb = options(1);
op_phi = options(2);
op_disp = options(3);
op_single = options(4);

if op_disp == 0 || op_disp == 1
    disp_flag = 0;
elseif op_disp == 2 && N ~= log2(size(p,1))
    disp_flag = 0;
else
    disp_flag = 1;
    if op_fb == 0
        fig_p = 200;
    else
        fig_p = 100;
    end
end

if N == 1
    [phi MIP0 MIP1 prob] = phi_comp(options,M,M,x1,p,b_table);
    Big_phi = phi;
    phi_m = phi*ones(1,3);
    prob_cell{1} = prob;
else

%% x0 bipartion
if op_fb == 0
    N_max = N-1; % not include the whole system
elseif op_fb == 1 || op_fb == 3
    N_max = N; % include the whole system
end
[x0_b1_o x0_b2_o N_b] = bipartition(M,N_max);

%% x1 data
C_x1 = cell(2^N-1,1);
k = 1;
for i=1: N
    C = nchoosek(M,i);
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j=1: N_C
        C_j = C(j,:);
        C_x1{k} = C_j;
        k = k + 1;
    end
end

MIP0 = cell(2^N-1,1); % MIP in x0
MIP1 = cell(2^N-1,1); % MIP in x1
phi = zeros(2^N-1,1); % small phis

prob = cell(2^N-1,1);
prob1 = cell(2^N-1,1);
prob2 = cell(2^N-1,1);

%% computing small phis
for i_C=1: 2^N-1
    C_j = C_x1{i_C}; % given data of x1
    if op_disp ~= 2
        fprintf('C=%s\n',mod_mat2str(C_j));
    end
    [phi(i_C) MIP0{i_C} MIP1{i_C} prob{i_C} prob1{i_C} prob2{i_C}] ...
        =  phi_comp(options,M,C_j,x1,p,b_table, x0_b1_o,x0_b2_o,N_b);
end

prob_cell{1} = prob;
prob_cell{2} = prob1;
prob_cell{3} = prob2;

phi_m = zeros(N,3); % cumulative sum

for i_C=1: 2^N-1
    C = C_x1{i_C};
    i = length(C);
    phi_m(i,1) = phi_m(i,1) + phi(i_C);
    phi_m(i,2) = phi_m(i,2) + phi(i_C)/nchoosek(N,i);
end

for i=1: N
    if i > 1
        phi_m(i,3) = phi_m(i-1,3) + phi_m(i,1);
    else
        phi_m(i,3) = phi_m(i,1);
    end
end

Big_phi = phi_m(end,3);

%% display
% index_vec = sort_index(N);
y_max = 0;
for i_C=1: 2^N-1
    y_max = max(y_max,max(prob{i_C}));
end
y_max = y_max + 0.02;

for i_C=1: 2^N-1
    C = C_x1{i_C};
    i = length(C);
    fprintf('C=%s: ',mod_mat2str(C));
    
    if i > 1 || op_single == 1
        x0_p1 = MIP0{i_C};
        x1_p1 = MIP1{i_C};
        
        sc = make_title(C,op_fb,M,x0_p1,x1_p1);
        fprintf('%s',sc);
    end
    
    if abs(phi(i_C)) < 10^-8
        phi(i_C) = 0;
    end
    fprintf(' phi%d=%f\n',i,phi(i_C));
end

if disp_flag == 1
    
    x = find(phi~=0);
    N_x = length(x);
    
    if op_fb == 3
        if 2^N < 9
            fig_max = 9;
        else
            fig_max = 16;
        end
        fig_vec = 1: fig_max;
    else
        if 2^N <= 8
            fig_max = 8;
            fig_co = 1; % 1 column
            fig_vec = 1: fig_max;
        else
            fig_max = 16;
            fig_co = 2; % 2 column
            fig_vec = [2*(1: fig_max/2)-1 2*(1:fig_max/2)];
        end
    end
    
    N_ir = 1;
    in_vec = 2^N-1: -1 :1 ;
    for i_C=2^N-1: -1 : 1
        C = C_x1{i_C};
        i = length(C);
        
        fig_i = floor((in_vec(i_C)-1)/fig_max);
        prob_all = prob{i_C};
        
        
        sC = make_title(C,op_fb);
        if op_fb == 3
            figure(fig_p+10+fig_i)
            sq = sqrt(fig_max);
            subplot(sq,sq,fig_vec(in_vec(i_C)-fig_i*fig_max)),imagesc(reshape(prob_all,[2^N 2^N]))
            title(sC);
        else
            figure(fig_p+10+fig_i)
            subplot(fig_max/fig_co,fig_co,fig_vec(in_vec(i_C)-fig_i*fig_max)),bar(prob_all)
            axis([-Inf Inf 0 y_max]);
            title(sC)
        end
        
        if i > 1 || op_single == 1
            x0_p1 = MIP0{i_C};
            x1_p1 = MIP1{i_C};
            if op_fb == 0 || op_fb == 3
                prob_prod = prob_prod_comp(prob1{i_C},prob2{i_C},M,x1_p1,op_fb);
            elseif op_fb == 1
                prob_prod = prob_prod_comp(prob1{i_C},prob2{i_C},M,x0_p1,op_fb);
            end
            
            sc = make_title(C,op_fb,M,x0_p1,x1_p1);
            if op_fb == 3
                figure(fig_p+20+fig_i)
                % fprintf('C=%s\n',mat2str(C));
                % fprintf('p1=%s\n',mat2str(prob1{i_C}));
                % fprintf('p2=%s\n',mat2str(prob2{i_C}));
                x = log2(length(prob_prod))/2;
                subplot(sq,sq,fig_vec(in_vec(i_C)-fig_i*fig_max)),imagesc(reshape(prob_prod,[2^x 2^x]))
                title(sc);
            else
                figure(fig_p+20+fig_i)
                subplot(fig_max/fig_co,fig_co,fig_vec(in_vec(i_C)-fig_i*fig_max)),bar(prob_prod)
                axis([-Inf Inf 0 y_max])
                title(sc)
            end
        end
        
        if phi(i_C) ~= 0
            sC = make_title(C,op_fb);
            sPhi = [sC,': \phi=',num2str(phi(i_C),4)];
            if op_fb == 3
                fig_pi = floor((N_ir-1)/fig_max);
                figure(fig_p+30+fig_pi)
                subplot(sq,sq,fig_vec(N_ir-fig_max*fig_pi)),imagesc(reshape(prob_all,[2^N 2^N]));
                title(sPhi);
            else
                fig_pi = floor((N_ir-1)/fig_max);
                figure(fig_p+30+fig_pi)
                subplot(fig_max/fig_co,fig_co,fig_vec(N_ir-fig_max*fig_pi)),bar(prob_all)
                axis([-Inf Inf 0 y_max]);
                title(sPhi)
            end
            if N_ir == N_x
                sy = ['\Phi=',num2str(Big_phi,4)];
                xlabel(sy)
            end
            N_ir = N_ir + 1;
        end
        
    end
    
end

for i=1: N
    fprintf('%d: phi_cum=%f phi_sum=%f phi_mean=%f\n',i,phi_m(i,3),phi_m(i,1),phi_m(i,2));
end

end

fprintf('Big phi=%f\n',Big_phi);

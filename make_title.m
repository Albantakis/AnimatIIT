function string = make_title(C,op_fb,M,x0,x1)

% op_fb = 0 (forward) ,1(backward), 3(both)

if nargin == 2
    if op_fb == 0
        string = ['F: [x0]=',mod_mat2str(C)];
    elseif op_fb == 1
        string = ['B: [x1]=',mod_mat2str(C)];
    elseif op_fb == 3
        string = ['FB: [x0]=',mod_mat2str(C)];
    end
else
    N = length(M);
    if op_fb == 0 || op_fb == 3
        x0_i = [];
        for i=1: length(C)
            x0_i = [x0_i x0(x0==C(i))];
        end
        x0_o = pick_rest(x0,x0_i);
        if isempty(x0_o) == 1
            s0 = mod_mat2str(x0_i);
        else
            s0 = mix_title(x0_i,x0_o);
        end
    elseif op_fb == 1
        s0 = mod_mat2str(x0);
    end
    s1 = mod_mat2str(x1);
    
    if length(x0) ~= N
        x0_r = pick_rest(M,x0);
        if op_fb == 0 || op_fb == 3
            x0_ri = [];
            for i=1: length(C)
                x0_ri = [x0_ri x0_r(x0_r==C(i))];
            end
            x0_ro = pick_rest(x0_r,x0_ri);
            if isempty(x0_ro) == 1
                s0_r = mod_mat2str(x0_ri);
            else
                s0_r = mix_title(x0_ri,x0_ro);
            end
        elseif op_fb == 1
            s0_r = mod_mat2str(x0_r);
        end
    else
        s0_r = '[]';
    end
    
    if length(x1) ~= N
        if op_fb == 0 || op_fb == 3
            x1_r = pick_rest(M,x1);
        elseif op_fb == 1
            x1_r = pick_rest(C,x1);
        end
        if isempty(x1_r) == 1
            s1_r = '[]';
        else
            s1_r = mod_mat2str(x1_r);
        end
    else
        s1_r = '[]';
    end
    
    if op_fb == 0
        string = ['F: [x0]-[x1]: ',s0,'-',s1,' & ',s0_r,'-',s1_r];
    elseif op_fb == 1
        string = ['B: [x0]-[x1]: ',s0,'-',s1,' & ',s0_r,'-',s1_r];
    elseif op_fb == 3
        string = ['FB: [x0]-[x1]: ',s0,'-',s1,' & ',s0_r,'-',s1_r];
    end
    
end

end

function s = mix_title(x_i,x_o)
s = '[';
for i=1: length(x_i)
    if i ~= length(x_i)
        s = [s int2str(x_i(i)) blanks(1)];
    else
        s = [s int2str(x_i(i))];
    end
end

s = [s blanks(1) '('];
for i=1: length(x_o)
    if i ~= length(x_o)
        s = [s int2str(x_o(i)) blanks(1)];
    else
        s = [s int2str(x_o(i)) ')]'];
    end
end
end

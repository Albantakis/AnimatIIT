function [string_p string] = make_title_fb(MIP,op_context,op_min)

if nargin < 3
    op_min = 0;
end

if op_min == 0
    xp_b1 = MIP{1,1};
    xp_b2 = MIP{2,1};
    x0_b1 = MIP{1,2};
    x0_b2 = MIP{2,2};
    xf_b1 = MIP{1,3};
    xf_b2 = MIP{2,3};
    x0 = combine(x0_b1,x0_b2);
    xp = combine(xp_b1,xp_b2);
    xf = combine(xf_b1,xf_b2);
    
    sp_1 = mod_mat2str(xp_b1);
    s0_1 = mod_mat2str(x0_b1);
    sf_1 = mod_mat2str(xf_b1);
    sp_2 = mod_mat2str(xp_b2);
    s0_2 = mod_mat2str(x0_b2);
    sf_2 = mod_mat2str(xf_b2);
    
    s0 = mod_mat2str(x0);
    sp = mod_mat2str(xp);
    sf = mod_mat2str(xf);
    
    if op_context == 0
        string{1} = [sp,'-',s0];
        string{2} = [s0,'-',sf];
        string{3} = [sp,'-',s0,'-',sf];
        string_p{1} = [sp_1,'-',s0_1,'x',sp_2,'-',s0_2];
        string_p{2} = [s0_1,'-',sf_1,'x',s0_2,'-',sf_2];
        string_p{3} = [sp_1,'-',s0_1,'-',sf_1,'x',sp_2,'-',s0_2,'-',sf_2];
    else
        string = [sp,'-',s0,'-',sf];
        string_p = [sp_1,'-',s0_1,'-',sf_1,'x',sp_2,'-',s0_2,'-',sf_2];
    end
else
    string = cell(3,1);
    string_p = cell(3,1);
    for bf = 1: 2
        xpf_b1 = MIP{1,1,bf};
        xpf_b2 = MIP{2,1,bf};
        x0_b1 = MIP{1,2,bf};
        x0_b2 = MIP{2,2,bf};
        x0 = combine(x0_b1,x0_b2);
        xpf = combine(xpf_b1,xpf_b2);
        
        s0_1 = mod_mat2str(x0_b1);
        s0_2 = mod_mat2str(x0_b2);
        s0 = mod_mat2str(x0);
        
        if bf == 1
            sp = mod_mat2str(xpf);
            sp_1 = mod_mat2str(xpf_b1);
            sp_2 = mod_mat2str(xpf_b2);
            string{bf} = [sp,' <- ',s0];
            string_p{bf} = ['(',sp_1,' <- ',s0_1,') x (',sp_2,' <- ',s0_2,')'];
        else
            sf = mod_mat2str(xpf);
            sf_1 = mod_mat2str(xpf_b1);
            sf_2 = mod_mat2str(xpf_b2);
            string{bf} = [s0,' -> ',sf];
            string_p{bf} = ['(',s0_1,' -> ',sf_1,') x (',s0_2,' -> ',sf_2,')'];
        end
    end
    string{3} = ['(',sp,'<-',s0,'->',sf,')'];
    string_p{3} = ['(',sp_1,'<-',s0_1,') x (',sp_2,'<-',s0_2,') & (',s0_1,'->',sf_1,') x (',s0_2,'->',sf_2,')'];
end

end
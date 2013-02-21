%% Old options:

op_single = network.options(4);     % just needed for console output
op_context = network.options(6);    % just needed in make_title
op_min = network.options(9);        % just needed in make_title
op_console = options(10);
op_big_phi = network.options(11);
op_normalize = options(13);         % big phi
op_normalize = network.options(14); % small phi

op_complex = network.options(15);

op_small_phi = network.options(16);
op_big_phi_dist = options(17);      
op_average = network.options(18);
op_parallel = in_options(19);


%% new options
op_parallel = in_options(1);
op_average = network.options(2);
op_complex = network.options(3);
op_small_phi = network.options(4);
op_big_phi = network.options(5);
op_normalize = network.options(6); % small phi
op_normalize = options(7);         % big phi
op_console = options(8);
op_parfor = network.options(9); % used by Animat program
op_strongconn = network.options(10);
op_extNodes = network.options(11);
op_single = network.options(12);    % just needed for console output == 1
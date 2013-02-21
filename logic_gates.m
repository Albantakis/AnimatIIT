function output = logic_gates(input,logic_type,in_noise)
% LOGIC_GATES the probability an element will turn on given the inputs
%
% OUTPUT = logic_gates(INPUT, LOGIC_TYPE)
%
% Given a binary input vector, input, and the logic type, logic_type, of
% the element of interest, this function will return the probability that
% the element of interest will be on

% 0 to .5
global noise;

if nargin == 3
    noise = in_noise;
end

% broken mechanism outputting 0 always
if logic_type == 0
    output = 0;
    
% AND
elseif logic_type == 1
    
    output = all(input);

% OR
elseif logic_type == 2

    output = any(input);

% XOR
elseif logic_type == 3

    output = sum(input) == 1;
    
% COPY    
elseif logic_type == 4

    output = input(1);

% NOT - TODO: check that element only has one input
elseif logic_type == 5
    
    output = ~input(1);

% NULL
elseif logic_type == 6

    output = .5;
    return;

% MAJORITY
elseif logic_type == 7  %Larissa: Majority should be more than half? 2 out of 4 is not majority, but is minority??
    
    N = length(input);
    output = (sum(input)/N > .5);

% MINORITY
elseif logic_type == 8
    
    N = length(input);
    output = (sum(input)/N < .5);
    
%PARITY
elseif logic_type == 9
    
    output = mod(sum(input),2) == 1;

%Linear Threshold unit, threshold is logic_type -10. I.e. '12' means threshold '2'    
elseif logic_type >= 10 && logic_type < 20
    
    output = sum(input) >= logic_type-10;

%below Threshold unit, threshold is logic_type -20. I.e. '22' means under threshold '2'    
elseif logic_type >= 20
    
    output = sum(input) < logic_type-20;   
        
end

output = abs(output - noise);
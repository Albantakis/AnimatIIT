% % Normalize by number of connections
% function Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2)
%  
%     Norm = length(numerator_part1)*length(denom_part2) + length(numerator_part2)*length(denom_part1);
% 
% end 
 
 
% Normalize by average mutual information between numerator and denominator elements (old way) 
function Norm = Normalization(xp_1,xp_2,numerator_part1,numerator_part2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(numerator_part1),length(xp_2)) + min(length(numerator_part2),length(xp_1));
else
    Norm = min(length(numerator_part1),length(xp_2)) + min(length(numerator_part2),length(xp_1)) ...
        + min(length(numerator_part1),length(xf_2)) + min(length(numerator_part2),length(xf_1));
end
end

% % Normalize by number of connections from numerator to denominator elements
% % (regardless of number of numerator elements)
% function Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2)
%  
%     Norm = ~isempty(numerator_part1)*numel(denom_part2) + ~isempty(numerator_part2)*numel(denom_part1);
% 
% end 

% function [ bytes ] = getMemSize( variable )
% props = properties(variable); 
% if size(props, 1) < 1, bytes = whos(isvarname(variable)); bytes = bytes.bytes;
% else %code of Dmitry
%   bytes = 0;
%   for ii=1:length(props)
%       currentProperty = getfield(variable, char(props(ii)));
%       s = whos(varname(currentProperty));
%       bytes = bytes + s.bytes;
%   end
% end
% end

function ss = getString_starts(p)

if p<=.001, ss = '***';
elseif p<=.01, ss = '**';
elseif p<=.05, ss = '*';
else, ss = ' n.s.';
end
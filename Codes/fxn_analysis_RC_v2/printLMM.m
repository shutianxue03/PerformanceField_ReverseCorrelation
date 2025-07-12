

function printLMM(data, varNames, model)
nvar = length(data);
switch nvar
    case 3, tbl = table(data{1}(:), data{2}(:), data{3}(:),'VariableNames',varNames);
    case 4, tbl = table(data{1}(:), data{2}(:), data{3}(:), data{4}(:),'VariableNames',varNames);
    case 5, tbl = table(data{1}(:), data{2}(:), data{3}(:), data{4}(:), data{5}(:),'VariableNames',varNames);
    case 6, tbl = table(data{1}(:), data{2}(:), data{3}(:), data{4}(:), data{5}(:), data{6}(:), 'VariableNames',varNames);
end

lme = fitlme(tbl, model);
slopes = lme.Coefficients.Estimate;
p = lme.Coefficients.pValue;

fprintf('Note that below is y=a+bx1+cx2\n')
fprintf('%s = %.4f%s + %.4f x %s%s + %.4f x %s%s\n', ...
    varNames{1}, slopes(1), getString_starts(p(1)), ...
    slopes(2), varNames{2}, getString_starts(p(2)), ...
    slopes(3),  varNames{3}, getString_starts(p(3)))


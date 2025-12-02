function [anova_text, tbl] = print_nANOVA(varNames, data, IVs, nsubj, lengMode)

% anova_text=print_nANOVA(varNames, data, IVs, nsubj, leng)
% varNames: cell
% data: a column vector
% IVs: a cell array, each vector is a column vector, same size as data
% nsubj: number of subj
% lengMode: text version; 
%    1=the concise version (for printng in the command window)
%    0=the simple version (for plotting in the figure)

%%
if nargin < 5, lengMode = 1; end % long version
nVars = length(varNames);

if isscalar(varNames)
    [~, tbl] = anovan(data, IVs, 'varnames', varNames, 'display', 'off');
else
    [~, tbl] = anovan(data, IVs,'model','interaction','varnames', varNames, 'display', 'off');
end

if nVars > 1
    nrows = nchoosek(nVars,2) + nVars;
else
    nrows = 1;
end

% sumsqr_error = tbl{nrows+2, 2};   % Sum of squares for error
sumsqr_total = tbl{nrows+3, 2};          % Total sum of squares

anova_text=[];
for aa=1:nrows
    varName = tbl{aa+1, 1};
    sumsqr = tbl{aa+1, 2};         % Sum of squares for Factor 1
    df1 = tbl{aa+1, 3};
    df2 = df1 * (nsubj-1);
    F = tbl{aa+1, 6};
    p = tbl{aa+1, 7};

    partial_eta2 = sumsqr/sumsqr_total;
    if lengMode % the concise version (for printng in the command window)
        anova_text=[anova_text, sprintf('%s: F(%d,%d)=%.2f, p=%.3f, partial-eta2=%.3f\n', varName, df1, df2, F,p, partial_eta2)];
    else % the simple version (for plotting in the figure)
        anova_text=[anova_text, sprintf('%s: F=%.2f, p=%.3f, eta^2_p=%.3f\n', varName, F,p, partial_eta2)];
    end
    %         anova_text=[anova_text, sprintf('%s%s\n', varName, getString_starts(p))];
end

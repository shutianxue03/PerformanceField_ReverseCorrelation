function cohenD = fxn_getES(x1, x2, type)
% fxn_getES  Effect size (Cohen's d)
%
% type:
%   'independent' : pooled SD (standard Cohen's d)
%   'paired'      : Cohen's dz = mean(diff)/std(diff)

if nargin < 3 || isempty(type)
    type = 'paired';
end

% ensure column vectors
x1 = x1(:);
x2 = x2(:);

switch lower(type)

    case 'independent'
        % remove NaNs independently (independent samples)
        x1 = x1(~isnan(x1));
        x2 = x2(~isnan(x2));

        n1 = numel(x1);
        n2 = numel(x2);
        if n1 < 2 || n2 < 2
            cohenD = NaN;
            return
        end

        v1 = var(x1, 0);  % sample variance
        v2 = var(x2, 0);

        sp2 = ((n1-1)*v1 + (n2-1)*v2) / (n1 + n2 - 2); % pooled variance
        if sp2 <= 0
            cohenD = NaN;
            return
        end

        cohenD = abs(mean(x1) - mean(x2)) / sqrt(sp2);

    case 'paired'
        % paired: keep only matched non-NaN pairs
        ok = ~isnan(x1) & ~isnan(x2);
        x1 = x1(ok);
        x2 = x2(ok);

        if numel(x1) < 2
            cohenD = NaN;
            return
        end

        d = x1 - x2;
        sd = std(d, 0);   % sample SD
        if sd == 0
            cohenD = NaN;
            return
        end

        cohenD = abs(mean(d)) / sd;  % Cohen's dz

    otherwise
        error('type must be ''independent'' or ''paired''.');
end
end
% ============================================================
% fxn_plotPF_diagram
% ------------------------------------------------------------
% Last updated by Shutian Xue on 07/20/2025

% Purpose:
%   Plot population field (PF) diagrams for different location groups
%   and marker sizes, saving the figures to specified folders.
%
% Inputs:
%   str_locgroup               - String label for the location group
%   sz_marker_all              - Array of marker sizes to plot
%   x_ref, y_ref               - Cell arrays of reference line coordinates
%   x_allLoc, y_allLoc         - Arrays of node coordinates
%   colors_allLoc              - Nx3 array of RGB colors for nodes
%   strParts_all               - Cell array of part names
%   indParts_all               - Cell array of indices for each part
%   outline_ecc                - Array of radii for outline arcs
%   sz_line                    - Line width for plot elements
%   nameFolder_Figures_diagrams- Output folder for saving figures

function fxn_plotPF_diagram(str_locgroup, sz_marker_all,  x_ref, y_ref, x_allLoc, y_allLoc, colors_allLoc, strParts_all, indParts_all, outline_ecc, sz_line, nameFolder_Figures_diagrams)
nParts = length(indParts_all);

for sz_marker = sz_marker_all

    if sz_marker == sz_marker_all(1), str_size = 'small'; else, str_size = 'big'; end
    for iPart = 1:nParts
        strParts = strParts_all{iPart};
        indParts = indParts_all{iPart};
        %%%%%%%%
        x_allLoc_part = x_allLoc(indParts);
        y_allLoc_part = y_allLoc(indParts);
        colors_allLoc_part = colors_allLoc(indParts, :);
        nNodes = length(x_allLoc_part);
        nRef = length(x_ref);

        figure('Position', [1e3 0 500 500], 'Color', 'w'); hold on;
        axis square, axis off;
        xlim([-2, 2])
        ylim([-2, 2])

        % Draw arcs
        theta = linspace(0, 2*pi, 200);
        for r = outline_ecc
            x = r * cos(theta);
            y = r * sin(theta);
            plot(x, y, 'k-', 'LineWidth', sz_line);
        end

        % Draw vertical dashed line
        for iRef = 1:nRef
            plot(x_ref{iRef}, y_ref{iRef}, 'k-', 'LineWidth', sz_line);
        end

        % Draw nodes
        for iNode = 1:nNodes
            x = x_allLoc_part(iNode);
            y = y_allLoc_part(iNode);
            plot(x, y, 'o', 'MarkerFaceColor', colors_allLoc_part(iNode, :), 'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker, 'LineWidth', sz_line)
        end

        % Define the name of the folder
        nameFolder = sprintf('%s/Size_%s', nameFolder_Figures_diagrams, str_size);
        if isempty(dir(nameFolder)), mkdir(nameFolder), end
        % Save the figure
        saveas(gcf, sprintf('%s/%s_%s.jpg', nameFolder, str_locgroup, strParts))

    end % iPart
    close all
end
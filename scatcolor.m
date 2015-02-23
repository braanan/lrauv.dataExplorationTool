function [cmap,coli,col] = scatcolor(x,range)
% Generate colormap for vector data
% x     = data vector
% range = color range (i.e 256)

% Some cool colormaps to play around with:
% linspecer(length(range))
% hsv(length(range))
% CT = cbrewer('div', 'Spectral', length(range));
cmap = colormap(jet(length(range))); % flipud(CT)

% colormapeditor
% pause

% color matrix for time series
coli = zeros(length(x),1);
col = zeros(length(x),3);
for j = 1:length(x)
    
    if ~isnan(x(j))
        if x(j)>=0    % case value is positive
            [~,coli(j)] = min(abs(abs(x(j))-range));
        elseif x(j)<0 % case value is negative
            [~,coli(j)] = min(abs(abs(x(j))+range));
        end
        
        col(j,:) = cmap(coli(j),:);
    
    else
        col(j,:) = NaN(1,3);
    end
end
end
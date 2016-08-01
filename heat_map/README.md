# This folder contains multiple block functions writing in Matlab to plot a [heat map] of a N*2 dataset

* [getBin](getBin.m)
* [calFreeEnergy](calFreeEnergy.m)
* [plot\_FreeEnergy_heatmap](plot_FreeEnergy_heatmap.m)

## Usage example:

```matlab
    Bins = getBin(data, BinWidth);
    F_energy = calFreeEnergy(Bins, T);
    [xGridNum, yGridNum] = size(Bins);
    x = data(:, 1);
    y = data(:, 2);
    X = linspace(min(x), max(x), xGridNum+1);
    Y = linspace(min(y), max(y), yGridNum+1);
    
    % Plot Free energy
    plot_FreeEnergy_heatmap(X, Y, Bins, F_energy);
    
    % Annotate Fig
    xlim = [min(x) max(x)];
    ylim = [min(y) max(y)];
    title(TitleName);
    set(gca,'YLim', ylim, 'XLim', xlim);
    set(gca, 'LineWidth', 2, 'FontSize', 14);
```

[heat map]:https://en.wikipedia.org/wiki/Heat_map

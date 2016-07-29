% plot the heatmap based on free energy
% X, Y are evenly spaced numbers over a specified interval

function plot_FreeEnergy_contour(X, Y, Bins, F_energy)

 % Plot contour
    [xGridNum, yGridNum] = size(Bins);
    for i=1:xGridNum
        for j=1:yGridNum
            if(Bins(i, j)>=1)  
                fill(X([i, i+1, i+1, i]), ...
                     Y([j, j, j+1, j+1]), ...
                     F_energy(i, j), ...
                     'linestyle', ...
                     'none');hold on
            end
        end
    end
    colormap jet
end 

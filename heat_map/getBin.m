%{ 
   Read in a 2D dataset contains X, Y data and the width of bin
   Return a Bin showing occurance of each grid
%}

function Bin = getBin(Data, xBinInfo, yBinInfo, xBinRange, yBinRange)
    x = Data(:, 1);
    y = Data(:, 2);
    if nargin == 3
        xRange = max(x) - min(x);
        yRange = max(y) - min(y);
        x_start = min(x);
        x_end = max(x);
        y_start = min(y);
        y_end = max(y);
    else
        x_start = xBinRange(1);
        x_end = xBinRange(2);
        y_start = yBinRange(1);
        y_end = yBinRange(2);
        xRange = x_end - x_start;
        yRange = y_end - y_start;
    end
    xisInt=~mod(xBinInfo, 1)
    yisInt=~mod(yBinInfo, 1)
    if xisInt | yisInt
        xGridNum = xBinInfo
        yGridNum = yBinInfo
    else
        xGridNum = int16(xRange/xBinInfo);
        yGridNum = int16(yRange/yBinWInfo);
    end
    Bin = zeros(xGridNum+1, yGridNum+1);
    [ConfNum, temp] = size(Data);
    for i=1:ConfNum
        % add 0.5 to convert starting index to 1 instead of 0
        xIndex = int16((Data(i, 1)-x_start)*xGridNum/xRange+0.5);
        yIndex = int16((Data(i, 2)-y_start)*yGridNum/yRange+0.5);
        Bin(xIndex, yIndex) = Bin(xIndex, yIndex) + 1;
    end
end

% Pis: absolute permutations
% Pijs: relative permutations
% I: edges
function [] = visualize_synth_match_data_all(Pis, Pijs, I, colors)

numCams = length(Pis);
N = size(Pis{1},2);
dimGrid = int32(sqrt(N));
dimGridBig = int32(sqrt(numCams));
e = length(I);

xys = {};

if (~exist('colors','var'))
    numColors = e;
    colors = distinguishable_colors(numColors);
end

[gx, gy] = meshgrid(1:dimGrid, 1:dimGrid);
gx = gx(:);
gy = gy(:);

cx = mean(gx);
cy = mean(gy);

coeff = 2;
[cx, cy] = meshgrid(coeff*cx*(1:dimGridBig), coeff*cy*(1:dimGridBig));
cx = cx(:);
cy = cy(:);

%clf;
hold on,
for k=1:numCams
    x = gx + cx(k);
    y = gy + cy(k);
    scatter(x,y);
    
    nums = (1:length(gx))';
    nums = num2str(nums);
    nums = cellstr(nums);
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    hold on,text(double(x)+dx, double(y)+dy, nums);
    
    xys = [xys [x y]];
end

axis equal;
colorId = 0;
for k=1:e
    i = I(1, k); j = I(2, k);
    
    if (i~=j) 
        
        Pi = Pis{i};
        Pj = Pis{j};
        Pij = Pijs{k};
        
        xyi = xys{i};
        xyj = xys{j};
        indxyiP = Pij * (1:length(xyi))';
        xyiP = xyj(indxyiP,:);
        
        hold on, line([xyi(:,1) xyiP(:,1)]', [xyi(:,2) xyiP(:,2)]', 'color', colors(mod(colorId, numColors)+1,:));
        colorId = colorId + 1;
    end
end

axis equal;

end
% Pis: absolute permutations
% Pijs: relative permutations
% points: points in the camera frame (will be shifted for better display)
% I: edges
% dims: x and y dimensions of the grid in points
function [] = visualize_synth_match_data(Pis, Pijs, points, I, dims)

n = length(Pis);
dimx = dims(1);
dimy = dims(2);
e = length(I);
offset = 0;

xys = {};

colors = distinguishable_colors(n);

%clf;
hold on,
for k=1:n
    pxy = points{k};
    x=double(pxy(:,1));
    y=double(pxy(:,2));
    x = offset+double(dimx*(k-1))*0.1+x;
    y = y - offset/3;
    scatter(x,y);
    
    nums = (1:length(pxy))';
    nums = num2str(nums);
    nums = cellstr(nums);
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    hold on,text(x+dx, y+dy, nums);
    offset = offset+3;
    
    xys = [xys [x y]];
end

axis equal;
colorId = 0;
for k=1:e
    i = I(1, k); j = I(2, k);
    
    if (i==j-1) % only plot the chain
        
        Pi = Pis{i};
        Pj = Pis{j};
        Pij = Pijs{k};
        
        xyi = xys{i};
        xyj = xys{j};
        indxyiP = Pij * (1:length(xyi))';
        xyiP = xyj(indxyiP,:);
        
        hold on, line([xyi(:,1) xyiP(:,1)]', [xyi(:,2) xyiP(:,2)]', 'color', colors(mod(colorId,n)+1,:));
        colorId = colorId + 1;
    end
end

axis equal;

end
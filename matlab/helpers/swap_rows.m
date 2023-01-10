function [P] = swap_rows(P, i, j)

rowi = P(i, :);
P(i,:) = P(j, :);
P(j, :) = rowi;

end
function [Ps] = perms_transform_first_to_eye(Ps)

K = length(Ps);

Ptrans = Ps{1}';
Ps{1} = eye(size(Ptrans,1));

for i=2:K
    Ps{i} = Ptrans*Ps{i};
end

end
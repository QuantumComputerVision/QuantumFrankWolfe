function [Pijs, I] = perms_Z_to_rel(Z, nViews)

k = nViews;
n = size(Z,1)./k;
Pijs={};
I=[];
for i=1:n
    for j=1:n
        Zsub = Z((i-1)*k+1:i*k, (j-1)*k+1:j*k);
        if(all(Zsub(:)==0))
            continue;
        end
        
        I=[I [i;j]];
        Pijs = [Pijs full(Zsub)];
    end
end

end
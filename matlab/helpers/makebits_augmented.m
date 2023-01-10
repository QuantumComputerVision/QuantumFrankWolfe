function A = makebits_augmented(n, dimIdentity)
    I = eye(dimIdentity);
    A = dec2bin(int32(2^n-1):int32(-1):int32(0))-'0'; 
    A = [repmat(I(:)', size(A,1), 1) A];
end
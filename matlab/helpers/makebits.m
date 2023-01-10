function A = makebits(n)
    A = dec2bin(2^n-1:-1:0)-'0'; 
end
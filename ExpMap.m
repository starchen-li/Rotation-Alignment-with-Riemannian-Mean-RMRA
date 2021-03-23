function P = ExpMap(B, S)
    BP = B^(1/2);
    BN = B^(-1/2);
    P  = BP * expm(BN * S * BN) * BP;
end

function S = LogMap(B, P)
    BP = B^(1/2);
    BN = B^(-1/2);
    S  = BP * logm(BN * P * BN) * BP;
end

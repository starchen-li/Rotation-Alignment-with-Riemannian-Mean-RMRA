function [CovsMT,D,M1,M2] =  Re-center(Covs)
M{1} = RiemannianMean(cat(3, Covs{:,1}));
M{2} = RiemannianMean(cat(3, Covs{:,2}));
D = RiemannianMean(cat(3, M{1}, M{2}));
M1 = M{1};
M2 = M{2};
M11 = M1^(-1/2);
M22 = M2^(-1/2);

for ss = 1 : 2
    for ii = 1 : size(Covs, 1)
        if ss == 1
            CovsMT{ii,ss} = M11*Covs{ii,ss}*M11;
        else
            CovsMT{ii,ss} = M22*Covs{ii,ss}*M22;
        end
       
    end
end
end

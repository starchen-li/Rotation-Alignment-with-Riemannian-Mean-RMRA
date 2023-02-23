function CovsUT =  RMRA(Covs,D0,M1,M2)
if nargin == 4
    D = D0;
    M{1} = M1;
    M{2} = M2;
else
    M{1} = RiemannianMean(cat(3, Covs{:,1}));
    M{2} = RiemannianMean(cat(3, Covs{:,2}));
    D = D0;
end
% D    = RiemannianMean(cat(3, M{1}, M{2}));
Mlabel_s = [1 2 3 4];
Mlabel_t = [1 2 3 4];

M_src =D;
M_tar1 = M{1};

options = struct('lambda',0.5,'dim',22,'kernel_type','primal','gamma',2,'data_feature',22); 
[M_src_new1,M_tar_new1,A1] = M-TCA(M_src,Mlabel_s',M_tar1,Mlabel_t',options);
B1 = A1^(1/2);
B1 = real(B1);
[U1,R] =qr(B1);

M_src = D;
M_tar2 = M{2};

options = struct('lambda',0.5,'dim',22,'kernel_type','primal','gamma',2,'data_feature',22);  
[M_src_new2,M_tar_new2,A2] = M-TCA(M_src,Mlabel_s',M_tar2,Mlabel_t',options);
B2 = A2^(1/2);
B2 = real(B2);
[U2,R] =qr(B2);


for ss = 1 : 2
    for ii = 1 : size(Covs, 1)
        if ss == 1
            U = U1;
            CovsUT{ii,ss} = U'*Covs{ii,ss}*U;
        else
            U = U2;
            CovsUT{ii,ss} = U'*Covs{ii,ss}*U;
        end
       
    end
end
end

function [X_src_new,X_tar_new,A] = M-TCA(X_src,Y_src,X_tar,Y_tar,options)
   %% Set options
    lambda = options.lambda;              
    dim = options.dim;                    
    kernel_type = options.kernel_type;    
    gamma = options.gamma; 
    data_feature = options.data_feature;
    
   %% Construct MMD matrix 
    X = [X_src',X_tar'];
    X = X*diag(sparse(1./sqrt(sum(X.^2))));
    [feature,n] = size(X);
    n = n/data_feature;
    ns = size(X_src,1)/data_feature;
    nt = size(X_tar,1)/data_feature;
    
    E = [];  %% Different from the original TCA
    for i = 1:ns
        E = [E;eye(feature)];
    end
	e = [1/ns*E;-1/nt*E];
% 	C = length(unique(Y_src));

	%%% M0
% 	M = e * e' * C;  %multiply C for better normalization
    M = e * e';
    
    M = M / norm(M,'fro');

	%% Centering matrix H, Different from the original TCA
    E = [E;E];
    a = eye(feature); 
    b = repmat({a},n,1);
    H = blkdiag(b{:});
    H = H - 1/n * (E*E');

	%% Calculation
    
	if strcmp(kernel_type,'primal')
        
	[A,~] = eigs(X*M*X'+lambda*eye(feature),X*H*X',dim,'SM');
        
    	Z = A'*X;
        Z = Z * diag(sparse(1./sqrt(sum(Z.^2))));
	X_src_new = Z(:,1:ns*feature)';
	X_tar_new = Z(:,ns*feature+1:end)';
    else
    	K = kernel_jda(kernel_type,X,[],gamma);
    	[A,~] = eigs(K*M*K'+lambda*eye(n*feature),K*H*K',dim,'SM');
    	Z = A'*K;
        Z = Z*diag(sparse(1./sqrt(sum(Z.^2))));
        X_src_new = Z(:,1:ns*feature)';
	X_tar_new = Z(:,ns*feature+1:end)';
    end

end


function K = kernel_jda(ker,X,X2,gamma)

    switch ker
        case 'linear'

            if isempty(X2)
                K = X'*X;
            else
                K = X'*X2;
            end

        case 'rbf'

            n1sq = sum(X.^2,1);
            n1 = size(X,2);

            if isempty(X2)
                D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X'*X;
            else
                n2sq = sum(X2.^2,1);
                n2 = size(X2,2);
                D = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X'*X2;
            end
            K = exp(-gamma*D); 

        case 'sam'

            if isempty(X2)
                D = X'*X;
            else
                D = X'*X2;
            end
            K = exp(-gamma*acos(D).^2);

        otherwise
            error(['Unsupported kernel ' ker])
    end
end

function miu = gmm2(X, K_or_centroids)
% ============================================================
% Expectation-Maximization iteration implementation of
% Gaussian Mixture Model.
%
% PX = GMM(X, K_OR_CENTROIDS)
% [PX MODEL] = GMM(X, K_OR_CENTROIDS)
%
%  - X: N-by-D data matrix.
%  - K_OR_CENTROIDS: either K indicating the number of
%       components or a K-by-D matrix indicating the
%       choosing of the initial K centroids.
%
%  - PX: N-by-K matrix indicating the probability of each
%       component generating each point.
%  - MODEL: a structure containing the parameters for a GMM:
%       MODEL.Miu: a K-by-D matrix.
%       MODEL.Sigma: a D-by-D-by-K matrix.
%       MODEL.Pi: a 1-by-K vector.
% ============================================================
 
    threshold = 1e-20;
    %数据是N个。每个数据有D维
    [N, D] = size(X);
 
    %表示k是一个标量还是一个矢量
    if isscalar(K_or_centroids)
        K = K_or_centroids;
        % randomly pick centroids
        rndp = randperm(N);
        centroids = X(rndp(1:K), :);
    else
        K = size(K_or_centroids, 1);
        centroids = K_or_centroids;
    end
 
    % initial values
    %对于EM算法
    %E就是计算后验，M就是是均值，方差最大化
    [pMiu pPi pSigma] = init_params();
 
    Lprev = -inf;
    while true
	%计算后验概率
        Px = calc_prob();
      
        % new value for pGamma 
	%repmat(pPi,N,1)相当于有N X 1个pPi
        pGamma = Px .* repmat(pPi, N, 1);
        %pGamma=r(nk)
	pGamma = pGamma ./ repmat(sum(pGamma, 2), 1, K);
 	
        % new value for parameters of each Component
        Nk = sum(pGamma, 1);
        pMiu = diag(1./Nk) * pGamma' * X;
        pPi = Nk/N;
	%update the gamma
        for kk = 1:K
            Xshift = X-repmat(pMiu(kk, :), N, 1);
            pSigma(:, :, kk) = (Xshift' * ...
                (diag(pGamma(:, kk)) * Xshift)) / Nk(kk);
        end
 
        % check for convergence
	%检测收敛
        L = sum(log(Px*pPi'));
        if L-Lprev < threshold
            break;
        end
        Lprev = L;
    end
 
%     nargout=0;
%     if nargout == 1
%         varargout = {Px};
%     else
        model = [];
        model.Miu = pMiu;
        model.Sigma = pSigma;
        model.Pi = pPi;
        varargout = {Px, model};
        miu=pMiu;
%     end
 
    %这个初始化不是太明白
    function [pMiu pPi pSigma] = init_params()
        pMiu = centroids;
        pPi = zeros(1, K);
        pSigma = zeros(D, D, K);
 
        % hard assign x to each centroids
        distmat = repmat(sum(X.*X, 2), 1, K) + ...
            repmat(sum(pMiu.*pMiu, 2)', N, 1) - ...
            2*X*pMiu';
	    %dummy是distmat的索引，labels是索引index,
        [dummy labels] = min(distmat, [], 2);
 
        for k=1:K
            Xk = X(labels == k, :);
            pPi(k) = size(Xk, 1)/N;
            pSigma(:, :, k) = cov(Xk);
        end
    end
 
    function Px = calc_prob()
    %有k个component!!
        Px = zeros(N, K);
        for k = 1:K
	    %remat pMiu(k, :)表示是一个gaussian model的矩阵
            Xshift = X-repmat(pMiu(k, :), N, 1);
	    %k个二维矩阵。pSigma。
            inv_pSigma = inv(pSigma(:, :, k));
	    %tmp=(x*inv_sigma*x)^2
            tmp = sum((Xshift*inv_pSigma) .* Xshift, 2);
            coef = (2*pi)^(-D/2) * sqrt(det(inv_pSigma));

	    %列向量 Px(:,K)=
            Px(:, k) = coef * exp(-0.5*tmp);
        end
    end
end

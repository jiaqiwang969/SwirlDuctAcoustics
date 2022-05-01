%% Clustering Using Gaussian Mixture Models
function [cutOffLine]=GMM_Cluster(lam,crLayer) %10 36 -2 35
mode= find(real(lam)>-50000&real(lam)<50000&imag(lam)>1&imag(lam)<5000&abs(lam).^2>1E-3);
X = [real(lam(mode)) imag(lam(mode))]; 
[n,p] = size(X);
rng(3); % For reproducibility

k = 2;
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
d = 500;
x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); % Increase number of EM iterations
gmfit = fitgmdist(X,k,'CovarianceType',Sigma{1},...
            'SharedCovariance',SharedCovariance{1},'Options',options);
gmfitX=[gmfit.mu(:,1)];
[qn1,qn2]=max([abs(crLayer(1)-gmfitX(1))+abs(crLayer(2)-gmfitX(1));abs(crLayer(1)-gmfitX(2))+abs(crLayer(2)-gmfitX(2))]);
cutOffLine=gmfitX(qn2);
%plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
end

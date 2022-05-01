%% Clustering Using Gaussian Mixture Models
%% How Gaussian Mixture Models Cluster Data
% Gaussian mixture models (GMM) are often used for data clustering.
% Usually, fitted GMMs cluster by assigning query data points to the
% multivariate normal components that maximize the component posterior
% probability given the data. That is, given a fitted GMM,
% <docid:stats_ug.brx2uny-1 cluster> assigns query data to the component
% yielding the highest posterior probability. This method of assigning a
% data point to exactly one cluster is called _hard_ clustering. For an
% example showing how to fit a GMM to data, cluster using the fitted model,
% and estimate component posterior probabilities, see
% <docid:stats_ug.bra9fvn Cluster Data from Mixture of Gaussian
% Distributions>.
%%
% However, GMM clustering is more flexible because you can view it as a
% _fuzzy_ or _soft clustering_ method. Soft clustering methods assign a
% score to a data point for each cluster. The value of the score indicates
% the association strength of the data point to the cluster. As opposed to
% hard clustering methods, soft clustering methods are flexible in that
% they can assign a data point to more than one cluster.  When clustering
% with GMMs, the score is the posterior probability.  For an example of
% soft clustering using GMM, see <docid:stats_ug.buqq83f Cluster Gaussian
% Mixture Data Using Soft Clustering>.
%%
% Moreover, GMM clustering can accommodate clusters that have different
% sizes and correlation structures within them. Because of this, GMM
% clustering can be more appropriate to use than, e.g, _k_-means
% clustering.
%%
% Like most clustering methods, you must specify the number of desired
% clusters before fitting the model. The number of clusters specifies the
% number of components in the GMM. For GMMs, it is best practice to also
% consider the: 
%
% * Component covariance structure.  You can specify diagonal or full
% covariance matrices, or whether all components have the same covariance
% matrix.
% * Initial conditions. The Expectation-Maximization (EM) algorithm fits
% the GMM. Like the _k_-means clustering algorithm, EM is sensitive to
% initial conditions and might converge to a local optimum.  You can
% specify your own starting values for the parameters, specify initial
% cluster assignments for data points or let them be randomly chosen, or
% specify to use the _k_-means ++ algorithm.
% * Regularization parameter. If, for example, you have more predictors
% than data points, then you can regularize for estimation stability.
%
%% Covariance Structure Options
% Load Fisher's iris data set.  Consider clustering the
% sepal measurements.
function GMM_Cluster_formal(lam) %10 36 -2 35
mode= find(real(lam)>-100&real(lam)<100&imag(lam)>-500&imag(lam)<500&real(lam).^2>1E-3);
X = [real(lam(mode)) imag(lam(mode))]; 
[n,p] = size(X);
rng(3); % For reproducibility

% figure;
% plot(X(:,1),X(:,2),'.','MarkerSize',15);
% title(num2str(m));
% xlabel('Sepal length (cm)');
% ylabel('Sepal width (cm)');
%% 
% The number of components, _k_, in a GMM determines number of
% subpopulations or clusters.  In this figure, it is difficult to determine
% whether two, three, or perhaps more components are appropriate. A GMM
% increases in complexity as _k_ increases.
%%
% Each component has a covariance matrix. Geometrically, the covariance
% structure determines the shape of a confidence ellipsoid drawn over a
% subpopulation or cluster.  You can specify whether the covariance
% matrices for all components are diagonal or full, or whether all
% components have the same covariance matrix.  Each combination of
% specifications determines the shape and orientation of the ellipsoids.
%%
% Fit GMMs to the data and examine the effects of specifying all
% combinations of covariance structure options on the shape of the
% ellipsoids. That is, specify all combinations of the name-value pair
% arguments |'CovarianceType'| and |'SharedCovariance'|. Covariance
% structure specifications apply to all components. For illustration,
% specify that there are three components.  To draw the ellipsoids:
% 
% # Use the fitted GMM to cluster a grid covering the plane composed of the
% extremes of the measurements.
% # Obtain the score that specifies a 99% probability threshold for each
% confidence region.  This specification determines the length of the major
% and minor axes of the ellipsoids.
% # Color the ellipse using a similar color to its cluster.
%
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

figure;
c = 1;
for i = 1:nSigma;
    for j = 1:nSC;
        gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i},...
            'SharedCovariance',SharedCovariance{j},'Options',options);
        clusterX = cluster(gmfit,X);
        mahalDist = mahal(gmfit,X0);
        subplot(2,2,c);
        h1 = gscatter(X(:,1),X(:,2),clusterX);
        hold on;
            for m = 1:k;
                idx = mahalDist(:,m)<=threshold;
                Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);

                uistack(h2,'bottom');
            end
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s, SharedCovariance = %s',...
            Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'1','2','3'});
        hold off
        c = c + 1;
    end
end
end

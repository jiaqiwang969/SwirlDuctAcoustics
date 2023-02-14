function [fitresult, gof] = createFit(r, u0)
%CREATEFIT(R,U0)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : r
%      Y Output: u0
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 11-Feb-2023 13:12:41 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( r, u0 );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );




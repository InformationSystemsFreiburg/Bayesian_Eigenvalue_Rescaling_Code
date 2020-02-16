function [CovBER,EigValBER] = BER(SampleRet,distribution,spacing,RollingCorr)
% Author: Torsten Mörstedt, Bernhard Lutz, Dirk Neumann @ University of
% Freiburg, Chair of Information Systems Research
% Function calculates the Bayesian Eigenvalue Rescaled Covariance as
% outlined in "Bayesian Eigenvalue Rescaling for Covariance Estimation in Portfolio Optimization"
% Setting used in the paper: [CovBER,EigValBER] = fnBER(ReturnMatrix,'Beta_2_6','Linear',HistoricRollingCorrelationVector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory Input:   
%       SampleRet: TxN matrix of asset returns
% Optional Input:
%       distribution: 'Normal','Beta_2_3','Beta_2_6'
%       spacing: 'Log','Linear','Wigner'
%       RollingCorr: 1xT matrix of rolling average cross-correlation values based on the
%       return matrix basket
% Example: 
%       [CovBER,EigValBER] = BER(randn(100,10),'Normal','Log',[0.4,0.35,0.3]);
%           randn(100,10):  TxN matrix of asset returns
%           [0.4,0.35,0.3]: 1xT matrix of rolling average cross-correlation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Management of input variables
    if nargin < 2
        distribution = 'Beta_2_6'; 
    end
    if nargin < 3
        spacing      = 'Linear';
    end
    if nargin < 4
        tempCorr    = corr(SampleRet);
        tempCorr(eye(size(tempCorr,1))) = NaN;
        RollingCorr = nanmean(nanmean(tempCorr));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BER Covariance
    CovSample                       = cov(SampleRet);
    [EigValBER,EigVecBER]           = BER_Inference(RollingCorr,CovSample,distribution,spacing); % Infere Bayesian eigenvalues  
    CovBER                          = EigVecBER*sqrt(diag(EigValBER))*inv(EigVecBER);            % Recreate (Bayesian Eigenvalue Rescaled) Covariance matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




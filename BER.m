function [CovBEF,EigValBEF] = BEF(SampleRet,distribution,spacing,RollingCorr)
% Author: Torsten Mörstedt, Bernhard Lutz, Dirk Neumann
% Function calculates the Bayesian Eigenvalue Rescaling Covariance Matrix as
% outlined in "Bayesian Bayesian Eigenvalue Rescaling for CovarianceEstimation in Portfolio Optimization"
% Setting used in the paper: [CovBEF,EigValBEF] = fnBEF(ReturnMatrix,'Beta_2_6','Linear',HistoricRollingCorrelationVector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mandatory Input:   
%       SampleRet: TxN matrix of asset returns
% Optional Input:
%       distribution: 'Normal','Beta_2_3','Beta_2_6'
%       spacing: 'Log','Linear','Wigner'
%       RollingCorr: 1xT matrix of rolling average cross-correlation values based on the
%       return matrix basket
% Example: 
%       [CovBEF,EigValBEF] = BEF(randn(100,10),'Normal','Log',[0.4,0.35,0.3]);
%           randn(100,10):  TxN matrix of asset returns
%           [0.4,0.35,0.3]: 1xT matrix of rolling average cross-correlation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Management of input variables
    if nargin < 2
        distribution = 'normal'; 
    end
    if nargin < 3
        spacing = 'Log';
    end
    if nargin < 4
        tempCorr = corr(SampleRet);
        tempCorr(eye(size(tempCorr,1))) = NaN;
        RollingCorr = nanmean(nanmean(tempCorr));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEF Covariance
    CovSample                       = cov(SampleRet);
    [EigValBEF,EigVecBEF]           = BEF_Inference(RollingCorr,CovSample,distribution,spacing); % Infere Bayesian eigenvalues  
    CovBEF                          = EigVecBEF*sqrt(diag(EigValBEF))*inv(EigVecBEF); % Recreate (Bayesian Eigenvalue Fitted) Covariance matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




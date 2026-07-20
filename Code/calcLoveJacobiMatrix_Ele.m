function J_Ele = calcLoveJacobiMatrix_Ele(freq,pvPred,DValue0,m0,DEN0,paraInd)
%   Summary of this function goes here.
%   J_Ele = calcLoveJacobiMatrix_Ele(freq,pvPred,DValue0,m0,DEN0,paraInd)
%   Detailed explanation goes here.
%   The function is for calculating the element of Jacobi matrix of the
%   misfit function of Love wave dispersion inversion.
%
%
%   IN
%           freq: the frequency vector, Hz.
%         pvPred: the phase velocity, m/s.
%        DValue0: the dispersion function value at the current model of the
%                 frequency freq.
%             m0: current model point, is consists of [S-wave velocity,
%                 layer thickeness].
%           DEN0: the density of current model point, kg/m^3.
%        paraInd: the index of model parameter of the Jacobi matrix element.
%
%
%  OUT
%          J_Ele: the element of Jacobi matrix.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2025-2030
%  Revision:  1.0  Date: 8/30/2025
%
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

delta = 1e-2;
if pvPred~=0
    mDelta = m0;
    mDelta(paraInd) = m0(paraInd)+delta;
    m = [mDelta,DEN0];
    DValue_right = fastlovecalc(pvPred,freq,m);
    D_m = (DValue_right-DValue0)/(delta);
    
    m = [m0,DEN0];
    DValue_right = fastlovecalc(pvPred+delta,freq,m);
    D_pv = (DValue_right-DValue0)/(delta);
    J_Element = -D_m/D_pv;
end

if pvPred==0
    J_Element = 0;
end

J_Ele = J_Element;
end


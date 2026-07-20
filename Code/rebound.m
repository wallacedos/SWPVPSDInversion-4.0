function mRb = rebound(m,Bounds)
%   Summary of this function goes here.
%   mRb = rebound(m,Bounds)
%   Detailed explanation goes here.
%   The function is for controlling the bounds of the input model 'm'.
%
%   IN
%              m: the input model, it is a row vector.
%         Bounds: the constraint matrix of the model bound, its size is n*2
%                 , n is the number of the parameters of model 'm'.
%
%
%  OUT
%            mRb: the output model after bound constraint, it is a row
%                 vector.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2024-2030
%  Revision:  1.0  Date: 2/26/2024
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

n = length(m);
n = length(m);
for i=1:n
    if m(i)<Bounds(i,1)
        m(i) = Bounds(i,1);
    end
    
    if m(i)>Bounds(i,2)
        m(i) = Bounds(i,2);
    end
end
mRb = m;
end


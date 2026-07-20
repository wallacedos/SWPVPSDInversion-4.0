function ccc = LoveSearchRoot(f,m,cc,n_mode)
%   Summary of this function goes here.
%   ccc = LoveSearchRoot(f,m,cc,n_mode)
%   Detailed explanation goes here.
%   The function is for searching the phase velocities of Love waves.
%
%   IN
%              f: the frequency, Hz.
%              m: the model vector, such as [S-wave velocity, 
%                 layer thickness.density].
%             cc: the trial phase velocity vector, it is a row vector,
%                 m/s.
%         n_mode: the mode number needs to compute.
%
%  OUT
%            ccc: the output phase velocity of Love waves, m/s.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2025-2030
%  Revision:  1.0  Date: 9/3/2025
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

nc = length(cc);
r = zeros(1,nc);
ccc = zeros(1,n_mode);
r(1)=fastlovecalc(cc(1),f,m);
n=1;
if(r(1)==0)
    ccc(n)=cc(1);
    n=n+1;
    if n==n_mode+1
        return;
    end
end
for j=2:nc
    r(j)=fastlovecalc(cc(j),f,m);
    if(r(j)==0)
        ccc(n)=cc(j);
        n=n+1;
    end
    if r(j)*r(j-1)<0
        ccd1=cc(j-1);ccd2=cc(j);rrd1=r(j-1);rrd2=r(j);
        for icd=1:8
            ccd=(ccd2+ccd1)/2;
            rrd=fastlovecalc(ccd,f,m);
            if rrd1*rrd<0
                ccd2=ccd;
                rrd2=rrd;
            end
            if rrd2*rrd<0
                ccd1=ccd;
                rrd1=rrd;
            end
        end
        ccc(n)=(ccd1+ccd2)/2;
        n=n+1;
    end
    if n==n_mode+1
        return;
    end
end
end


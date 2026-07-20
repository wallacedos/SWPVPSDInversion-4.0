function y= fastlovecalc(x,f,m)
%   Summary of this function goes here
%   function y=fastlovecalc(x,f,m)
%   The function is the key function for calculating the disperive function
%   of Love wave.
%   Detailed explanation goes here
%   The function is called by the function named alclovemulti.
%
%IN
%    x: The trail phase velocity, m/s.
%    f: frequency, Hz.
%    m: the model vector, it is a row vector, such as [S-wave velocity,
%       layer thickness, density]. 
%
%   References: Schwab F A,Knopoff L.Fast surface wave and free mode
%               computions[M].In:Methods in Computational Physics.Bolt B
%               A(ed.).New York: Academic Press,1972:87-180.
%
%       History
%          1) The function name has changed from fastlovecalc(x) to 
%             fastlovecalc(x, f, m) 9/3/2025.
%
%  Author(s): Yan Yingwei
%  Copyright: 2017-2030 
%  Revision: 2.0  Date: 9/3/2025.
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

n = length(m);
n = floor((n+1)/3+1e-6);
VS = m(1:n);
H = m(n+1:2*n-1);
den = m(2*n:3*n-1);
k=2*pi*f/x;
s=-den(n)*VS(n)^2*sqrt((1-x^2/VS(n)^2));
T0=[s -1i];
F=T0;
a=zeros(2,2,n-1);
for i=1:n-1
    r=-1i*sqrt(1-x^2/VS(i)^2);
    h=H(i);
    u=den(i)*VS(i)^2;
    Q=k*r*h;
    a(:,:,i)=[cos(Q) 1i*sin(Q)/(u*r);1i*u*r*sin(Q) cos(Q)];
end
for i=n-1:-1:2
    F=F*a(:,:,i);
end
i=1;
a1=zeros(2,1);
a1(1,1)=a(1,1,i);
a1(2,1)=a(2,1,i);
F=F*a1;
y=real(F);
end


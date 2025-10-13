function p=v_vonmisespdf(x, mu, k)
%0002 %V_VONMISESPDF Von Mises probability distribution P=(x,m,k)
%0003 %
%0004 %  Inputs:  X         matrix of input values (in radians)
%0005 %           M         mean angle of distribution (in radians)
%0006 %           K         concentration parameter
%0007 %
%0008 % Outputs:  P         matrix of probability density values (same size as X)
%0009 %                     (with no output argument, the function will plot a graph)
%0010 %
%0011 % The von Mises distribution describes the pdf of an angle over the range [0,2pi).
%0012 % For large K, the distribution approximates a Gaussian with mean M and
%0013 % variance 1/K. For small K, the distribution is uniform.
%0014
%0015 %      Copyright (C) Mike Brookes 1997-2011
%0016 %      Version: $Id: v_vonmisespdf.m 10865 2018-09-21 17:22:45Z dmb $
%0017 %
%0018 %   VOICEBOX is a MATLAB toolbox for speech processing.
%0019 %   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%0020 %
%0021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0022 %   This program is free software; you can redistribute it and/or modify
%0023 %   it under the terms of the GNU General Public License as published by
%0024 %   the Free Software Foundation; either version 2 of the License, or
%0025 %   (at your option) any later version.
%0026 %
%0027 %   This program is distributed in the hope that it will be useful,
%0028 %   but WITHOUT ANY WARRANTY; without even the implied warranty of
%0029 %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%0030 %   GNU General Public License for more details.
%0031 %
%0032 %   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>0
    p=exp(k*cos(x-mu))/(2*pi*besseli(0,k));
else
    if nargin<1 || isempty(x)
        x=linspace(-pi,pi,181);
    end
    if nargin<2 || isempty(mu)
        mu=0;
    end
    if nargin<3 || isempty(k)
        k=[0 pow2(-1:3)];
    end
    nm=length(mu);
    nk=length(k);
    np=max(nm,nk);
    pp=zeros(length(x),np);
    pl=cell(np,1);
    for i=1:np
        mi=mu(1+rem(i-1,nm));
        ki=k(1+rem(i-1,nk));
        pp(:,i)=v_vonmisespdf(x(:),mi,ki);
        pl{i}=sprintf('\\mu=%.1f, \\kappa=%.1f',mi,ki);
    end
    plot(x,pp);
    v_axisenlarge([-1 -1.05]);
    legend(pl,'location','northeast');
end
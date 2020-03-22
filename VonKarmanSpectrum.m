function [S] = VonKarmanSpectrum(f,U,stdVel,L,component)
% ---------------------------------------------
% INPUT
% f: float; frequency is [1 x Nfreq]
% V: float; Mean wind speed Normal to the deck is [1 x 1]
% stdVel : float; std of speed is [1 x 1]
% L=  float; turbulence length scales is [1 x 1]
% Iturb = ; float; turbulence intensity is [1 x 1]
% component : string; is 'u','v' or 'w'
% ---------------------------------------------
% OUTPUT
% S: float; [Nzz x Nyy] matrix 
% ---------------------------------------------
% ---------------------------------------------
% Author: E. Cheynet 02.11.2016
% 
%%

% Von Karman coefficent
a=[-5/6, -11/6];
% dimension of output
%calculation of S/std^2
f_hat = L.*U.^(-1).*f;

if strcmpi(component,'u'),
    S =  U.^(-1).*4.*L.*stdVel.^2.*(1+70.7.*f_hat.^2).^(a(1));
elseif or(strcmpi(component,'w'),strcmpi(component,'v')),
    S=  U.^(-1).*4.*L.*stdVel.^2.*(1+70.7.*4.*f_hat.^2).^(a(2)).*(1+188.4.*4*f_hat.^2);
else
    error('component unknown')
end













end


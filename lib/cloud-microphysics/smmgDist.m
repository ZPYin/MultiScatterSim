function [n_r] = smmgDist(r, varargin)
% SMMGDIST single-mode modified gamma distribution for cloud droplet size spectrum.
%
% USAGE:
%    [n_r] = smmgDist(r)
%
% INPUTS:
%    r: numeric
%        radius. (¦Ìm)
%
% KEYWORDS:
%    kappa: numeric
%        shape parameter. (default: 7)
%    sigma_e: numeric
%        extinction coefficient of cloud. (default: 2e-2 m^{-1})
%    R_eff: numeric
%        effective radius. (¦Ìm)
%
% OUTPUTS:
%    n_r: numeric
%        cloud droplet number concentration at given radius. (m^{-3}¦Ìm^{-1})
%
% REFERENCES:
%    [1] Donovan D, Klein Baltink H, Henzing J, De Roode S, Siebesma A. A
%        depolarisation lidar-based method for the determination of liquid-cloud
%        microphysical properties. Atmospheric Measurement Techniques. 2015;8(1)
%        :237--66.
%    [2] Ahmad W, Zhang K, Tong Y, Xiao D, Wu L, Liu D. Water Cloud Detection 
%        with Circular Polarization Lidar: A Semianalytic Monte Carlo Simulation
%        Approach. Sensors. 2022;22(4):1679.
%
% HISTORY:
%    2023-05-22: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'r', @isnumeric);
addParameter(p, 'kappa', 7, @isnumeric);
addParameter(p, 'sigma_e', 6e-2, @isnumeric);
addParameter(p, 'R_eff', 9, @isnumeric);

parse(p, r, varargin{:});

R_m = p.Results.R_eff ./ (p.Results.kappa + 2);   % Eq.(8) Ref (1)
k = (p.Results.kappa + 1) .* p.Results.kappa ./ (p.Results.kappa + 2).^2;   % Eq.(11) Ref (1)
N_d = 1 / (2 * pi) * p.Results.sigma_e ./ (k .* p.Results.R_eff.^2);   % Eq.(4) Ref(2)

n_r = N_d ./ R_m ./ factorial(p.Results.kappa - 1) .* ...
      (r ./ R_m).^(p.Results.kappa - 1) .* exp(- r ./ R_m) * 1e12;

end
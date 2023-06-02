function [S1, S2, Qext, Qsca] = MieCalc(m, x, theta, varargin)
% MIECALC Mie Calculation.
%
% USAGE:
%    [S1, S2, Qext, Qsca] = MieCalc(m, x, theta)
%
% INPUTS:
%    m: complex number
%        refractive index.
%    x: numeric
%        size parameter.
%    theta: numeric
%        scattering angle. [degree]
%
% OUTPUTS:
%    S1: complex number
%        S1 parameter. (see Ref (1))
%    S2: complex number
%        S2 parameter. (see Ref (1))
%    Qext: numeric
%        extinction efficiency.
%    Qsca: numeric
%        scattering efficiency.
%
% REFERENCES:
%    [1] Ahmad W, Zhang K, Tong Y, Xiao D, Wu L, Liu D. Water Cloud Detection
%        with Circular Polarization Lidar: A Semianalytic Monte Carlo Simulation
%        Approach. Sensors. 2022;22(4):1679.
%
% HISTORY:
%    2023-05-23: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'm', @isnumeric);
addRequired(p, 'x', @isnumeric);
addRequired(p, 'theta', @isnumeric);
addParameter(p, 'debug', false, @islogical);

parse(p, m, x, theta, varargin{:});

if (x == 0)                 % To avoid a singularity at x=0

    S1 = 0;
    S2 = 0;
    Qext = 0;
    Qsca = 0;

elseif (x > 0)              % This is the normal situation

    nmax = round(2 + x + 4 * x.^(1/3));
    n1 = nmax - 1;
    n = (1:nmax);
    cn = 2 * n + 1;
    c1n = n .* (n + 2) ./ (n + 1);
    c2n = cn ./ n ./ (n + 1);
    x2 = x .* x;

    f = Mie_ab(m, x);
    anp = (real(f(1,:))); 
    anpp = (imag(f(1,:)));
    bnp = (real(f(2,:)));
    bnpp = (imag(f(2,:)));
    dn = cn .* (anp + bnp);

    q = sum(dn);
    Qext = 2 * q / x2;

    en = cn .* (anp .* anp + anpp .* anpp + bnp .* bnp + bnpp .* bnpp);
    q = sum(en);
    Qsca = 2 * q / x2;

    an = f(1, :);
    bn = f(2, :);
    pt = Mie_pt(cos(theta / 180 * pi), nmax);
    pin = pt(1, :);
    tin = pt(2, :);
    n2 = (2 * n + 1) ./ (n .* (n + 1));
    pin = n2 .* pin;
    tin = n2 .* tin;
    S1 = (an * pin' + bn * tin');
    S2 = (an * tin' + bn * pin');

end
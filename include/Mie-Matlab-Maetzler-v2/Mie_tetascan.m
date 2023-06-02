function result = Mie_tetascan(m, x, nsteps, varargin)
% Computation and plot of Mie Power Scattering function for given 

% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 

% according to Bohren and Huffman (1983) BEWI:TDD122

% C. Mätzler, May 2002.

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'm', @isnumeric);
addRequired(p, 'x', @isnumeric);
addRequired(p, 'nsteps', @isnumeric);
addParameter(p, 'visible', 'off', @ischar);

parse(p, m, x, nsteps, varargin{:});

m1=real(m); m2=imag(m);

nx=(1:nsteps); dteta=pi/(nsteps-1);

teta=(nx-1).*dteta;

for j = 1:nsteps

    u=cos(teta(j));

    a(:,j)=Mie_S12(m,x,u);

    SL(j)= real(a(1,j)'*a(1,j));

    SR(j)= real(a(2,j)'*a(2,j));

end

y=[teta teta+pi;SL SR(nsteps:-1:1)]'; 

if strcmpi(p.Results.visible, 'on')
    polar(y(:,1),y(:,2));

    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x));

    xlabel('Scattering Angle');
end

result=y; 


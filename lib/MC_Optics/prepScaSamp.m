function optics = prepScaSamp(optics, mixture, tailDegrees)
%PREPSCASAMP P11 proposal with exact discrete p/q weights.
% The mixture only changes sampling variance, never the physical phase matrix.
if nargin<2,mixture=.2;end
if nargin<3,tailDegrees=170;end
assert(isscalar(mixture) && mixture>=0 && mixture<1);
assert(isscalar(tailDegrees) && tailDegrees>=0 && tailDegrees<=180);
if isfield(optics,'proposalCorrection')
    assert(isequal(size(optics.proposalCorrection),size(optics.cdfS11)));
    assert(all(isfinite(optics.proposalCorrection),'all') && all(optics.proposalCorrection>=0,'all'));
    assert(all(optics.cdfS12==0,'all'), 'External scalar proposals require a zero CDF12.');
else
    mass=diff([zeros(size(optics.cdfS11,1),1),optics.cdfS11],1,2);
    assert(all(isfinite(mass),'all') && all(mass>=-1e-12,'all'));
    mass=max(mass,0);
    total=sum(mass,2);assert(all(total>0));
    tail=double(reshape(optics.theta,1,[])>=tailDegrees);
    probability=sum(mass.*tail,2)./total;
    effectiveMix=mixture*double(probability>0);
    q=(1-effectiveMix).*mass + effectiveMix.*(mass.*tail)./max(probability,realmin);
    assert(max(abs(sum(q,2)./total-1))<1e-10);
    optics.proposalCorrection=mass./max(q,realmin);
    optics.cdfS11=cumsum(q,2);optics.cdfS12(:)=0;
end
end

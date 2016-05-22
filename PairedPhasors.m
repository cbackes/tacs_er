function out = PairedPhasors(th,rho)
% function c

% assumes samples > nCategories.

if ~all(size(th)==size(rho))
    error('sizes of inputs do not match')
end

if size(th,1)<size(th,2)
    th = th'; rho = rho';
end
N           = size(th,1);
nCategories = size(th,2);

out = [];
out.Categ = [];
out.Combs = [];
out.nCategories = nCategories;
nCombs = nchoosek(nCategories,2);
out.nCombs = nCombs;
combFields  = {'CatCombs','DeltaMTh','DeltaRho','DeltaMZTh','DeltaRayPZ'};
catFields   = {'meanTh','meanRho','meanZTh','RayleighPZ'};

% means by category
meanTh      = zeros(nCategories,2);
meanRho     = zeros(nCategories,1);
meanZTh     = zeros(nCategories,2);
RayleighPZ  = zeros(nCategories,2);
for jj = 1:nCategories
    thj     = th(:,jj);
    rhoj    = rho(:,jj);
    
    % angle nums
    meanTh(jj,1)  = angle(mean(exp(1j*thj)));
    meanTh(jj,2)  = abs(mean(exp(1j*thj)));
    
    % radix nums
    meanRho(jj)     = mean(rhoj);
    
    % combined
    z = mean(rhoj.*exp(1j*thj));
    meanZTh(jj,1)       = angle(z);
    meanZTh(jj,2)       = abs(z);
    
    % Rayleigh Test
    [RayleighPZ(jj,1),RayleighPZ(jj,2)] = circ_rtest(thj);
end

% save into output struct
for kk = 1:numel(catFields);
    out.Categ.(catFields{kk}) = eval(catFields{kk});
end

%%
% get across differences between categories.
CatCombs    = zeros(nCombs,2);
DeltaMTh    = zeros(nCombs,2);
DeltaRho    = zeros(nCombs,1);
DeltaMZTh   = zeros(nCombs,2);
DeltaRayPZ  = zeros(nCombs,2);

cnt= 1;
for jj = 1:nCategories
    thj     = th(:,jj);
    rhoj    = rho(:,jj);
    zj      = rhoj.*exp(1j*thj);
    for ii = (jj+1):nCategories
        thi     = th(:,ii);
        rhoi    = rho(:,ii);
        zi      = rhoi.*exp(1j*thi);
        CatCombs(cnt,:)     = [jj,ii];
        
        % angle nums
        DeltaMTh(cnt,1)     = angle(mean(exp(1j*(thj-thi))));
        DeltaMTh(cnt,2)     = abs(mean(exp(1j*(thj-thi))));
        
        % rho nums
        DeltaRho(cnt)       = mean(rhoj-rhoi);
        
        % combined
        DeltaMZTh(cnt,1)    = angle(mean(zi-zj));
        DeltaMZTh(cnt,2)    = abs(mean(zi-zj));
        
        % RayleighTest
        [DeltaRayPZ(cnt,1), DeltaRayPZ(cnt,2) ] = circ_rtest(thj-thi);
        cnt = cnt+1;
    end
end

% save into output struct
for kk = 1:numel(combFields);
    out.Combs.(combFields{kk}) = eval(combFields{kk});
end


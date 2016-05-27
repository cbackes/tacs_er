function simulations2(N)
%%
N=1000;

% Number of samples
nn = 300;

% Number of offset phases
No = 10;

% Visual Presentation Error Term
er = 2.5/(2*pi);

% von misses distribution
vm = @(k,th)(exp(k*cos(th))/(2*pi*besseli(0,k)));

th = 0:2*pi/5:2*pi-0.1;

offsetSet = linspace(0,2*pi,No);

da = zeros(N,No);
for ii = 1:N
    for jj = 1:No
        offset  = offsetSet(jj);
        th1     = th-offset;
        th2     = th-dth-offset;
        
        rho1    = vm(kp,th1);
        rho2    = vm(kp,th2);
        
        th1_r = randsample(th,nn,true,rho1);
        th1_r = th1_r(:)+unifrnd(-er,er,[nn,1]);
        th2_r = randsample(th,nn,true,rho2);
        th2_r = th2_r(:)+unifrnd(-er,er,[nn,1]);
        
        a1 = angle(mean(exp(1j*th1_r)));
        a2 = angle(mean(exp(1j*th2_r)));
        da(ii,jj) = mod(a1-a2,2*pi);
    end
end

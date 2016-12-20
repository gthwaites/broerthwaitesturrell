
% Main program in the new formulation
% Created 15/04/16 by Arthur Turrell
% Based on an earlier code by Greg Thwaites

% Hi both, just confirming that git is working by adding an inane comment.
% I will upload this to the server.

% In this code, variables per capita and absolute variables have the same
% units as the number of households is made up of a continuum which is
% normalised to one.
clear all
% close all

%% Computational variables
ninterp=2*10^2;           % number of points of wealth distribution at which policy functions evaluated
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter',1500,'MaxFunEvals',200);

%% Global constants
alpha=1/3;
bettadisp=.3;
bettamean=0.5
bettavec=linspace(bettamean-bettadisp,bettamean+bettadisp,3);                       % Discount rate on future actions

%% Solve tree model with constant beta
betta=bettamean;
% 1=(1-beta)*(p+(1-alpha))
% p=(1-betta)^-1-(1-alpha);
% 1=(1-beta)*(1+p)
p=betta/(1-betta);
r=alpha/p;
R=1+r;

%% Solve tree model with random beta, no frictions
bmax=(alpha+p)*5;                                       % Highest considered bequest 
abvec=linspace(0,bmax,ninterp);
mu=0;
kappa=0;

% check that function solves
[kexcess,a_ss,ab_ss,abmapmat2,k_ss]=tree(mu,kappa,abvec,alpha,bettavec,p,R);

fun = @(x) tree_clearing(mu,kappa,abvec,alpha,bettavec,x,R)
[x,fval,exitflag] = fmincon(fun,[p],[],[],[],[],[0],[Inf]);
        
% plot capital discrepancy as a function of R excess
parray=linspace(0,bmax,20);
pexcesstest=zeros(size(parray));
for ii=1:length(parray)
    pexcesstest(ii)=tree_clearing(mu,kappa,abvec,alpha,bettavec,parray(ii),R);
end
plot(parray,pexcesstest,parray,zeros(size(pexcesstest)))


% solve model for different alpha and betadisp, calc inequality
bettadispvec=linspace(0,.25,3);
alphavec=linspace(0.1,0.5,3);
lorenz=zeros(length(abvec),length(bettadispvec),length(alphavec));
gini=zeros(length(bettadispvec),length(alphavec));

for ii=1:length(bettadispvec)
    for jj=1:length(alphavec)
        bettadisp=bettadispvec(ii);
        alpha=alphavec(jj);
        bettavec=linspace(bettamean-bettadisp,bettamean+bettadisp,3);                       % Discount rate on future actions
        fun = @(x) tree_clearing(mu,kappa,abvec,alpha,bettavec,x,R)
        [x,fval,exitflag] = fmincon(fun,[p],[],[],[],[],[0],[Inf]);
        [kexcess,a_ss,ab_ss,abmapmat2,k_ss]=tree(mu,kappa,abvec,alpha,bettavec,p,R);
        
        totalwealth(ii,jj)=abvec*ab_ss;
        lorenz(:,ii,jj)=cumsum(ab_ss.*abvec')/totalwealth(ii,jj);
        gini(ii,jj)=(.5-trapz(cumsum(ab_ss),lorenz(:,ii,jj)))/.5;

    end
end


% compare wealth distribution for different values of alpha
%% Solve tree model with fixed cost and financial friction
kappa=0.1;
mu=0.2;

% plot capital discrepancy as a function of R excess
parray=linspace(x*.01,10*x,10);
Rarray=linspace(R*.01,3*R,20);
kexcesstest=zeros(length(parray),length(Rarray));
aexcesstest=zeros(length(parray),length(Rarray));
for ii=1:length(parray)
    for jj=1:length(Rarray)
        [kexcesstest(ii,jj),aexcesstest(ii,jj)]=tree(mu,kappa,abvec,alpha,bettavec,parray(ii),Rarray(jj));
    end
end
test=kexcesstest.^2+aexcesstest.^2;
[C,I]=min(test);
[C2,I2]=min(C);
pstart=parray(I(I2))
Rstart=Rarray(I2)
%pstart=10^-1
%Rstart=10^-1
fun = @(x) tree_clearing_kr(mu,kappa,abvec,alpha,bettavec,x(1),x(2))
[x,fval,exitflag]= fmincon(fun,[pstart Rstart],[],[],[],[],[0 0],[Inf Inf]);
        
return


%% Solve productive capital model for Rexcess given r
kappaarray=linspace(0.01,0.2,3);
muarray=linspace(0.1,0.25,5);
xarray=zeros(length(kappaarray),length(muarray));
totalwealth=xarray;
gini=xarray;
lorenz=zeros(length(abvec),length(kappaarray),length(muarray));
Rstart=.1;

% plot capital discrepancy as a function of R excess
Rexcessarray=linspace(-r,20,100);
for ii=1:length(Rexcessarray)
    [kexcesstest(ii),~,~,~,ktest(ii),k_sstest(ii)] = mktclearing_tree(r,Rexcessarray(ii),abvec,bettavec,0.1,0.1,alpha,delta);
end

[hAx,hLine1,hLine2] = plotyy(Rexcessarray,ktest,Rexcessarray, k_sstest)

title('Capital demand and supply')
xlabel('excess return on capital')

ylabel(hAx(1),'capital implied by excess return') % left y-axis
ylabel(hAx(2),'capital supplied by savers') % right y-axis

% return

for kk=1:length(kappaarray)
    for jj=1:length(muarray)
        
        fun = @(x) t_clearing(r,x,abvec,bettavec,kappaarray(kk),muarray(jj),alpha,delta);
        [xarray(kk,jj),fvalarray(kk,jj),exitflagarray(kk,jj)] = fmincon(fun,[2],[],[],[],[],[0],[Inf]);
        Rstart=xarray(kk,jj);
        [kexcess,a_ss,ab_ss] = mktclearing_tree(r,xarray(kk,jj),abvec,bettavec,kappaarray(kk),muarray(jj),alpha,delta);
        
        totalwealth(kk,jj)=abvec*ab_ss;
        lorenz(:,kk,jj)=cumsum(ab_ss.*abvec')/totalwealth(kk,jj);
        gini(kk,jj)=(.5-trapz(cumsum(ab_ss),lorenz(:,kk,jj)))/.5;
        R=1+xarray(kk,jj)+r;
        rho=(R-(1+r)*(1-muarray(jj)))/muarray(jj);
        k=((R+delta-1)/alpha)^(1/(alpha-1));
        w=(1-alpha)*k^alpha;
        cutoffwealth(:,kk,jj)=(kappaarray(kk)/w)./(1-((1+r)/rho).^bettavec)-1; % scale by wages
        
        % what %ile of the wealth distribution is a capitalist?
        for aa=1:size(bettavec)
            capshare(aa)=(abvec>=cutoffwealth(aa,kk,jj))*ab_ss;
            
        end
        capsharemat(kk,jj)=mean(capshare);
        
        % 90/10
    end
end
figure
subplot(2,2,1)
surf(muarray,kappaarray,xarray)
xlabel('\mu')
ylabel('\kappa')
zlabel('Excess return on capital')

subplot(2,2,2)
surf(muarray,kappaarray,gini)
xlabel('\mu')
ylabel('\kappa')
zlabel('Gini coefficient of wealth')

subplot(2,2,3)
%surf(muarray,kappaarray,log10(sqrt(fvalarray)))
surf(muarray,kappaarray,(fvalarray))
xlabel('\mu')
ylabel('\kappa')
zlabel('Capital discrepancy at optimum')

subplot(2,2,4)
surf(muarray,kappaarray,squeeze(mean(cutoffwealth,1)))
xlabel('\mu')
ylabel('\kappa')
zlabel('Mean cutoff value of wealth')

[~,~,ab_ss,abmapmat2] = mktclearing_tree(r,xarray(2,2),abvec,bettavec,kappaarray(2),muarray(2),alpha,delta);
        
figure
surf(abmapmat2)
save matlab.mat
return
%% Solve productive capital model for Rexcess with endogenous r
kappaarray=linspace(0.01,0.2,5);
muarray=linspace(0.1,0.5,10);
xarrayr=zeros(2,length(kappaarray),length(muarray));
totalwealthr=zeros(length(kappaarray),length(muarray));
ginir=totalwealthr;
lorenzr=zeros(length(abvec),length(kappaarray),length(muarray));
rstart=0;
Rexcessstart=.01;
% use previous result as new starting value
for kk=1:length(kappaarray)
    for jj=1:length(muarray)
        
        % loop over starting values r from -1 to 1 Rexcess from 0.01 to 1
        % evaluate and choose the lowest value
        
        fun = @(x) rk_clearing(x(1),x(2),abvec,bettavec,kappaarray(kk),muarray(jj),alpha,delta);
        [xarrayr(:,kk,jj),fvalarrayr(kk,jj),exitflagarrayr(kk,jj)]= fmincon(fun,[rstart Rexcessstart],[],[],[],[],[-1 0],[Inf Inf]);
        rstart=xarrayr(1,kk,jj);
        Rexcessstart=xarrayr(2,kk,jj);
        [kexcess,a_ss,ab_ss] = mktclearing_k_scalar(xarrayr(1,kk,jj),xarrayr(2,kk,jj),abvec,bettavec,kappaarray(kk),muarray(jj),alpha,delta);
        
        totalwealthr(kk,jj)=abvec*ab_ss;
        lorenzr(:,kk,jj)=cumsum(ab_ss.*abvec')/totalwealthr(kk,jj);
        ginir(kk,jj)=(.5-trapz(cumsum(ab_ss),lorenzr(:,kk,jj)))/.5;
        
        R=1+xarrayr(1,kk,jj)+xarrayr(2,kk,jj);
        r=xarrayr(1,kk,jj);
        rho=(R-(1+r)*(1-muarray(jj)))/muarray(jj);
        k=((R+delta-1)/alpha)^(1/(alpha-1));
        w=(1-alpha)*k^alpha;
        cutoffwealthr(:,kk,jj)=(kappaarray(kk)/w)./(1-((1+r)/rho).^bettavec)-1; % scale by wages
        
        % what %ile of the wealth distribution is a capitalist?
        for aa=1:size(bettavec)
            capsharer(aa)=(abvec>=cutoffwealthr(aa,kk,jj))*ab_ss;            
        end
        capsharematr(kk,jj)=mean(capsharer);
        
    end
end
figure
subplot(2,2,1)
surf(muarray,kappaarray,squeeze(xarrayr(1,:,:)))
xlabel('\mu')
ylabel('\kappa')
zlabel('Return on bonds, end r')

figure
subplot(2,2,2)
surf(muarray,kappaarray,squeeze(xarrayr(2,:,:)))
xlabel('\mu')
ylabel('\kappa')
zlabel('Excess return on capital, end r')

subplot(2,2,3)
surf(muarray,kappaarray,(fvalarrayr))
xlabel('\mu')
ylabel('\kappa')
zlabel('Capital discrepancy at optimum')

figure
surf(muarray,kappaarray,ginir)
xlabel('\mu')
ylabel('\kappa')
zlabel('Gini coefficient of wealth, end r')

return
%% Solve housing model for excess return and rent given r
totalwealthh=zeros(length(kappaarray),length(muarray));
ginih=totalwealthh;
xarrayh=zeros(2,length(kappaarray),length(muarray));
lorenzh=zeros(length(abvec),length(kappaarray),length(muarray));
w=1;
phi=0.2;
for kk=1:length(kappaarray)
    for jj=1:length(muarray)
        
        fun = @(x) Hh_clearing(r,x(1),x(2),abvec,bettavec,kappaarray(kk),muarray(jj),phi);
        [xarrayh(:,kk,jj)] = fmincon(fun,[.1 .1],[],[],[],[],[0 0],[Inf Inf]);
        [Hexcess,hexcess,a_ss,ab_ss] = mktclearing_h_scalar(r,xarrayh(1,kk,jj),xarrayh(2,kk,jj),abvec,bettavec,kappaarray(kk),muarray(jj),phi);
        
        totalwealthh(kk,jj)=abvec*ab_ss;
        lorenzh(:,kk,jj)=cumsum(ab_ss.*abvec')/totalwealthh(kk,jj);
        ginih(kk,jj)=(.5-trapz(cumsum(ab_ss),lorenzh(:,kk,jj)))/.5;
        
        % 90/10
    end
end
figure
surf(muarray,kappaarray,squeeze(xarrayh(1,:,:)))
xlabel('\mu')
ylabel('\kappa')
zlabel('Excess return on housing')

figure
surf(muarray,kappaarray,ginih)
xlabel('\mu')
ylabel('\kappa')
zlabel('Gini coefficient of housing wealth')

save results.mat
%% Now solve housing model for r, excess return and rent

% wealth distribution depends on cutoff and excess return

% other salient statistics
% fraction of households owning capital



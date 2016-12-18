function [F] = tree_clearing_kr(mu,kappa,abvec,alpha,bettavec,p,R)

[kexcess,a_ss,~] = tree(mu,kappa,abvec,alpha,bettavec,p,R);

F=kexcess^2+a_ss^2;

end

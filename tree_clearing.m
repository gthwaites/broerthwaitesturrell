function [F] = tree_clearing(mu,kappa,abvec,alpha,bettavec,p,r)

[kexcess,~,~] = tree(mu,kappa,abvec,alpha,bettavec,p,r);

F=kexcess;

end
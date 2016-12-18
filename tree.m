function [kexcess,a_ss,ab_ss,abmapmat2,k_ss]=tree(mu,kappa,abvec,alpha,bettavec,p,r)
%   Detailed explanation goes here

        abmapmat=zeros(length(abvec),length(abvec),length(bettavec));
        abhatmat=zeros(length(abvec),length(bettavec));
        abvecsquare=repmat(abvec,[length(abvec) 1]); % different starting wealth in the columns
        abvecprimesquare=repmat(abvec',[1 length(abvec)]); % different end wealth in the rows

        w=(1-alpha); % wage from this capital stock
        sigma=1;
        
if mu==0
    if kappa==0
            %% no borrowing, no fixed costs
        for kk=1:length(bettavec)
            r=alpha/p;
            ca=abvecsquare*(1+r)+w-abvecprimesquare;
            
            if sigma==1
                ua=(1-bettavec(kk))*log(abvecprimesquare)+bettavec(kk)*log(max(0,ca));
            else
                ua=(1-bettavec(kk))*((abvecprimesquare).^(1-sigma))/(1-sigma)+bettavec(kk)*(max(0,ca).^(1-sigma))/(1-sigma);
            end

            [Ma,Ia]=max(ua,[],1);
            decision=full(sparse(Ia,1:length(abvec),ones(length(abvec),1)));
            abmapmat(1:size(decision,1),1:size(decision,2),kk)=decision; % check that first two indices are right way round
            kmat(:,kk)=abvec(Ia);  %k_k=(betta/mu)*(abvec+w-kappa);
            amat(:,kk)=0;
            abmapmat2=mean(abmapmat,3); % average over betta to get the average transition

            % should we transpose this matrix or not?
            % Find fixed point of the wealth distribution - construct Markov transition matrix by integrating over abhatmat
            [eigenvectors,eigenvalues]=eig(abmapmat2);
            eigenvaluesr=round(10^5*eigenvalues)/10^5;
            candidates=eigenvectors(:,find(diag(eigenvaluesr)==1));
            candidates=candidates./repmat(sum(candidates,1),[size(candidates,1) 1]);
            k_cand=zeros(size(candidates,2),1);
            for jj=1:size(candidates,2)
                k_cand(jj)=mean(candidates(:,jj)'*kmat); % demand for trees
            end
            [kexcess,I]=min(abs(k_cand-1));
            a_ss=0; % should equal zero in equilibrium
            k_ss=k_cand(I);
            kexcess=k_ss-1;
            ab_ss=candidates(:,I);

            % check disconnect between kmat and abmat

            % calculate steady state capital for different steady
            % states, choose wealth distribution closest to starting
            % value
        end
    else % kappa>0, mu=0
            % does this case make sense? If kappa>0, then very poor
            % households might choose to lend. But if there is no limit to
            % borrowing then they can all gear up infinitely.
            % So no, it doesn't make sense
       % display an error
    end
    
else % mu>0 
    if kappa==0
        % everyone gears up to the same extent, no heterogeneity so no borrowing
        % display an error
    else
        rho=(1+alpha/p-(1+r)*(1-mu))/mu; % geared excess return on capital
        for kk=1:length(bettavec)
            
            ca=abvecsquare*(1+r)+w-abvecprimesquare;
            cb=abvecsquare*(1+rho)+w-kappa*p-abvecprimesquare;
            if sigma==1
                ua=(1-bettavec(kk))*log(abvecprimesquare)+bettavec(kk)*log(max(0,ca));
                ub=(1-bettavec(kk))*log(abvecprimesquare)+bettavec(kk)*log(max(0,cb));
            else
                ua=(1-bettavec(kk))*((abvecprimesquare).^(1-sigma))/(1-sigma)+bettavec(kk)*(max(0,ca).^(1-sigma))/(1-sigma);
                ub=(1-bettavec(kk))*((abvecprimesquare).^(1-sigma))/(1-sigma)+bettavec(kk)*(max(0,cb).^(1-sigma))/(1-sigma);
            end

            [Ma,Ia]=max(ua,[],1);
            [Mb,Ib]=max(ub,[],1);
            bchoose=(Ma<Mb);    % 1 if that point in abvec is a capitalist

            decision=full(sparse(Ia.*(1-bchoose)+Ib.*bchoose,1:length(abvec),ones(length(abvec),1)));
            abmapmat(1:size(decision,1),1:size(decision,2),kk)=decision; % check that first two indices are right way round
            kmat(:,kk)=bchoose.*abvec(Ib);  %k_k=(betta/mu)*(abvec+w-kappa);
            amat(:,kk)=(1-bchoose').*abvec(Ia)'/(1+r)-(1-mu)*kmat(:,kk);
        
            abmapmat2=mean(abmapmat,3); % average over betta to get the average transition

            % should we transpose this matrix or not?
            % Find fixed point of the wealth distribution - construct Markov transition matrix by integrating over abhatmat
            [eigenvectors,eigenvalues]=eig(abmapmat2);
            eigenvaluesr=round(10^5*eigenvalues)/10^5;
            candidates=eigenvectors(:,find(diag(eigenvaluesr)==1));
            candidates=candidates./repmat(sum(candidates,1),[size(candidates,1) 1]);
            k_cand=zeros(size(candidates,2),1);
            for jj=1:size(candidates,2)
                k_cand(jj)=mean(candidates(:,jj)'*kmat); % demand for trees
            end
            [kexcess,I]=min(abs(k_cand-1));
            k_ss=k_cand(I);
            kexcess=k_ss-1;
            
            a_ss=0; % should equal zero in equilibrium
            
            ab_ss=candidates(:,I);

            % check disconnect between kmat and abmat

            % calculate steady state capital for different steady
            % states, choose wealth distribution closest to starting
            % value
 end
end

end





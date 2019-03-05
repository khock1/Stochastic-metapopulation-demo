function [populations] = pop_dyn( patches, timesteps, populations, connct, dstbmort, natgrowth, natsettlmort, natadultmort, K, compt, extinctthresh  )

%run population dynamics

for t=2:timesteps
    
    %determine number of juveniles produced
    juv=zeros(patches,1);
    for p=1:patches
        juv(p,1)=populations(p,t-1)*natgrowth(p);
    end
    
    %use connectivity to determine number of available juveniles at a patch
    tmpJ=zeros(patches,patches);
    shutdown=[timesteps/2 timesteps timesteps];
    for p1=1:patches
        for p2=1:patches
            if t<shutdown(p1)
                tmpJ(p1,p2)=(connct(p1,p2,t)*juv(p1,1));
            else
                if p1~=p2
                    tmpJ(p1,p2)=(connct(p1,p2,t)*juv(p1,1));
                end
            end
        end
    end
    tmpS=zeros(patches,1);
    for p=1:patches
        tmpS(p,1)=sum(tmpJ(:,p));
    end

    %implement mortality
    tmpR=zeros(patches,1);
    for p=1:patches
        gamma=(1-natsettlmort(p,t));
        beta=(1/K(p,1))*((gamma/natadultmort(p,1))-(1/natgrowth(p)));
        tmpR(p,1)=(gamma*tmpS(p,1))/(1+(beta*tmpS(p,1)));
        if tmpR(p,1)<extinctthresh%population crashes completely if less than this number
            tmpR(p,1)=0;
        end
    end


    tmpP=tmpR+populations(:,t-1);
    
    %work out the popualtion size in next tiem step from birth, immigration, and mortality
    for p=1:patches
        if tmpP(p)>K(p,1)
            tmpP(p)=K(p1,1);
        end
        tmpP(p)=tmpP(p)*(1-natadultmort(p));
        tmpP(p)=tmpP(p)*(1-dstbmort(p,t));

    end
    populations(:,t)=tmpP;
    
end

end
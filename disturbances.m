function [dstbcmlmort, natsettlmort, dstbdistb] = disturbances(patches, timesteps, dodstb, climateC, natsettlmort, anthroWQ, pcycl, ...
    cyccatW, cycextnW, pwet, pblch, blchsevrW, blchextnW, pcots, cotsyrlmort, indvprobs, fixedcycl, fixedblch, fixedcots, fixedwety)

%Run different kinds fo disturbances to get cumulative mortality


%Cyclones---------------------------------------------------------------
ccatmortscale=[0.045 0.205 0.515 0.74 0.95];%mortality from categories
patchcyc=zeros(3,timesteps);
if indvprobs==0
    cyccat=zeros(1,timesteps);
else
    cyccat=zeros(patches,timesteps);
end
cycextn=zeros(1,timesteps);
cycyr=zeros(1,timesteps);
if dodstb(1)==1%determine whether cyclone occurs and at what intensity
    if indvprobs==0
        if fixedcycl==0
            for t=1:timesteps
                if rand<pcycl
                    cycyr(1,t)=1;
                    cyccat(1,t)=datasample(1:5,1,'weights',cyccatW);
                    cycextn(1,t)=datasample(1:4,1,'weights',cycextnW);
                end
            end
        else
            cycyr(1,fixedcycl(:,1))=1;
            cyccat(1,fixedcycl(:,1))=fixedcycl(:,2);
            cycextn(1,fixedcycl(:,1))=fixedcycl(:,3);
        end
    else
        for t=1:timesteps
            for p=1:patches
                if rand<pcycl(p,t)
                    cycyr(1,t)=1;
                    cyccat(p,t)=cyccatW(p,t);
                end
            end
        end
    end
end


%Floods and human impacts from floods on WQ--------------------------------
precip=zeros(1,timesteps);
cycrain=0.15;
settlmorts=[(1-0.441) (1-0.691) (1-0.941)];
wetyr=zeros(1,timesteps);
if dodstb(4)==1
    for t=1:timesteps
        if climateC(4)==1
            pwet=[(1/3-0.15) 1/3 (1/3+0.15)];%reduces chance of a wet year
        end
        tpwet=pwet;
        if cyccat(1,t)>0%if cycloen occurs, flood occurs
            tpwet(1)=pwet(1)+cycrain;
            tpwet(3)=pwet(3)-cycrain;
        end
        if anthroWQ==1%human impacts on WQ
            wsm=[100 50 1];
            wqeff=wsm+(t*0.5);
            tsettlmorts=0.941-(0.005*wqeff);
        else
            tsettlmorts=settlmorts;
        end
        if fixedwety==0
            precip(1,t)=datasample(1:3,1,'weights',tpwet);
        else
            precip(1,fixedwety(:,1))=fixedwety(:,2);
        end
        if precip(1,t)==1
            wetyr(1,t)=1;
            natsettlmort(3,t)=tsettlmorts(1);
        elseif precip(1,t)==2
            natsettlmort(3,t)=tsettlmorts(2);
        else
            wetyr(1,t)=-1;
            natsettlmort(3,t)=tsettlmorts(3);
        end
    end
end


%Bleaching------------------------------------------------------------------
blchmortscale=[(1-0.5066) (1-0.2567) (1-0.1300) (1-0.0659)];%from Hughes for tabular Acropora, for 2 4 6 8 weeks
patchblch=zeros(3,timesteps);
blyr=zeros(1,timesteps);
dryeffect=[-.15 0.1 0.05];%changes in bleaching extent prob due to dry year, more likely patch affected if no clouds
blchsevr=zeros(1,timesteps);
blchextn=zeros(1,timesteps);
if dodstb(2)==1
    if length(pblch)==1
        if climateC(2)==1
            pblch=0.1:pblch:0.55;
            pblch=repelem(pblch,round(timesteps/length(pblch)));
        else
            pblch(1,1:timesteps)=0.1;
        end
    end
    for t=1:timesteps
        if fixedblch==0
            if precip(1,t)==3%if it is a dry year, more likely bleaching will cover more patches
                tblchextnW=blchextnW+dryeffect;
            else
                tblchextnW=blchextnW;
            end
            if rand<pblch(t)
                blyr(1,t)=1;
                blchsevr(1,t)=datasample(1:4,1,'weights',blchsevrW);
                blchextn(1,t)=datasample(1:3,1,'weights',tblchextnW);
            end
        else
            blyr(1,fixedblch(:,1))=1;
            blchsevr(1,fixedblch(:,1))=fixedblch(:,2);
            blchextn(1,fixedblch(:,1))=fixedblch(:,3);
        end
        if blyr(t)==1
            if blchextn(1,t)==1
                patchblch(2,t)=blchmortscale(blchsevr(1,t));
            elseif blchextn(1,t)==2
                patchblch(2:3,t)=blchmortscale(blchsevr(1,t));
            else
                patchblch(:,t)=blchmortscale(blchsevr(1,t));
            end
        end
    end
end


%Calcualte mortalities----------------------------------------------------
if dodstb(1)==1
    for t=1:timesteps
        if blyr(t)~=1
            if cycyr(1,t)==1
                if indvprobs==0
                    if cycextn(1,t)==1
                        patchcyc(2,t)=ccatmortscale(cyccat(1,t));
                    elseif cycextn(1,t)==2
                        patchcyc(1:2,t)=ccatmortscale(cyccat(1,t));
                    elseif cycextn(1,t)==3
                        patchcyc(2:3,t)=ccatmortscale(cyccat(1,t));
                    else
                        patchcyc(:,t)=ccatmortscale(cyccat(1,t));
                    end
                else
                    for p=1:patches
                        if cyccat(p,t)>0
                            patchcyc(p,t)=ccatmortscale(cyccat(p,t));
                        end
                    end
                end
            end
        else
            cycyr(1,t)=0;
        end
    end
end


%COTS outbreaks------------------------------------------------------------
patchcots=zeros(3,timesteps);%not linked to coral cover
cotsyr=zeros(1,timesteps);
cotstimer=0;
if dodstb(3)==1
    for t=1:timesteps
        if fixedcots==0
            if cotstimer==0%possible to get outbreak immediately after the previous one has finished, not linked to coral abundance
                if precip(1,t)==3
                    tpcots=pcots*1.5;
                else
                    tpcots=pcots;
                end
                if rand<tpcots
                    cotstimer=1;%COTS outbreaks last 6 years
                    cotsyr(1,t)=1;
                end
            end
        elseif ismember(t,fixedcots)
            cotsyr(1,t)=1;
            cotstimer=1;
        end
        if cotstimer>0
            patchcots(2,t)=cotsyrlmort;
            cotstimer=cotstimer+1;
            if cotstimer>6
                cotstimer=0;
            end
        end
    end
end


%Add up mortalities from all disturbances----------------------------------
dstbcmlmort=zeros(3,timesteps);
for t=1:timesteps
    for p=1:size(patchcots,1)
        cml1=(1-patchcots(p,t));
        if patchcyc(p,t)>0
            cml1=cml1*(1-patchcyc(p,t));
        end
        if patchblch(p,t)>0
            cml1=cml1*(1-patchblch(p,t));
        end
        dstbcmlmort(p,t)=1-cml1;
    end
end

dstbdistb=vertcat(cycyr,blyr,cotsyr,wetyr);

end

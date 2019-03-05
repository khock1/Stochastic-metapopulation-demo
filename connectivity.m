function nets3 = connectivity(  type, patches, timesteps, biases )

%determine connectivity pattern

if biases==0
    biases=[0.7 0.1 0.01;
            0.01 0.7 0.1;
            0.01 0.1 0.7];
end

switch type
    case 'sr'%simple adjacency matrix
        nets3=repmat([1 0 0;0 1 0;0 0 1], [1 1 timesteps]);
    case 'fx'%fixed connectivity
        nets3=repmat(biases, [1 1 timesteps]);
    case 'vr'%uniformly variable connectivity using biases as a base
        nets3=zeros(patches,patches, timesteps);
        for t=1:timesteps
            for p=1:patches
                r = rand(1, patches);
                r=r.*biases(p,:);
                r = r / sum(r);
                nets3(p,1:patches,t)=r(1,1:patches);
            end
        end
    case 'wn'%add white noise to the connectvity pattern
        snr=0.2;%lower means more noise
        noise4selfrecr=1;
        nets3=zeros(patches,patches, timesteps);
        for p1=1:patches
            for p2=1:patches
                if noise4selfrecr==1
                    nets3(p1,p2,:)=awgn(repmat(biases(p1,p2),[timesteps 1]),snr,'measured');
                else
                    if p1==p2
                        nets3(p1,p2,:)=repmat(biases(p1,p2),[timesteps 1]);
                    else
                        nets3(p1,p2,:)=awgn(repmat(biases(p1,p2),[timesteps 1]),snr,'measured');
                    end
                end
            end
        end
        nets3(nets3<0)=0;
    case 'iwn'%add white noise to each link independently
        snr=[10 2 2];%lower means more noise
        noise4selfrecr=[0 0 0];
        limit=500;
        nets3=zeros(patches,patches, timesteps);
        for p1=1:patches
            for p2=1:patches
                if noise4selfrecr(p1)==1
                    if p1==p2
                        nets3(p1,p2,:)=awgn(repmat(biases(p1,p2),[timesteps 1]),snr(p1),'measured');
                    else
                        nets3(p1,p2,1:(timesteps-limit))=awgn(repmat(biases(p1,p2),[timesteps-limit 1]),snr(p1),'measured');
                    end
                else
                    if p1==p2
                        nets3(p1,p2,:)=repmat(biases(p1,p2),[timesteps 1]);
                    else
                        nets3(p1,p2,1:(timesteps-limit))=awgn(repmat(biases(p1,p2),[timesteps-limit 1]),snr(p1),'measured');
                    end
                end
            end
        end
        nets3(nets3<0)=0;
end
end

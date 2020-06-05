function Times=FiringTimes(l)
% computes the origin firing times given the firing propensities
[forks oris]=size(l);
for fork=1:forks
    U=rand(1,oris);
    for ori=1:oris
        if l(fork,ori)~=0
            Times(fork,ori)=-log(U(ori))/(l(fork,ori));
        else
            Times(fork,ori)=inf;
        end
    end
end

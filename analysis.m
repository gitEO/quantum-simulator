clear;
figure
fig;
test=0;
limiter=0.01;
n=1;
m=0;
for i=1:size(scale,2)
    if EnergyProb(i)<=limiter
        if test==1;
            test=0;
            a(n)=i;
            n=n+1;
        end
    else
        if test==0;
            test =1;
            a(n)=i;
            n=n+1;
            m=m+1;
            prob(m)=0;
        end
        prob(m)=prob(m)+EnergyProb(i);
    end
end



totProb=0;
for k=1:(size(a,2)/2)
    range=a(2*k-1):a(2*k);
    medio=0;
    under=0;
    totProb=totProb+prob(k);
    for j=range
        medio=medio+scale(j)*EnergyProb(j);
        under=under+EnergyProb(j);
    end
    energy(k)=medio/under;
end
energy
prob
totProb

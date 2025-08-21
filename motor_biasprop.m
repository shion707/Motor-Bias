function l=motor_biasVindras2005(data1,x)
% Prop vec

Ttheta=[0:pi/12:2*pi-0.001];
[tx,ty] = pol2cart(Ttheta,1);
startp=[x(1);x(2)];
T=[tx; ty];
planp=T-startp;
[ttemp,rtemp] = cart2pol(planp(1,:),planp(2,:));

f=find(ttemp<-0);
ttemp(f)=ttemp(f)+2*pi;

rtemp=x(3)*rtemp;
ttemp=ttemp+x(4)+rtemp;

error_theta=ttemp-Ttheta;

SE=(error_theta-data1);
SE=SE(:);
Sigm=nanstd(SE);
LLH=log(normpdf( SE,0 ,Sigm));
l=-nansum(nansum(LLH));






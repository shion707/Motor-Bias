function l=motor_biasVindras2005_VB(data1,x)
% Prop Vec + TG model

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

visualbias=[0:x(5) x(5)-0.01:-2*x(5)/(91-2*x(5)-1):-x(5) -x(5):-0.1]*x(6);
visualbias=[visualbias visualbias visualbias visualbias];

vb=visualbias(1:15:360);

error_theta=error_theta+vb;
f=find(error_theta<-pi);
error_theta(f)=error_theta(f)+2*pi;

SE=(error_theta-data1);
SE=SE(:);
Sigm=nanstd(SE);
LLH=log(normpdf( SE,0 ,Sigm));
l=-nansum(nansum(LLH));



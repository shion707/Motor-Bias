function l=motor_biasallp_vb(x,data1,data2)
% TG model

visualbias=[0:x(6) x(6)-0.01:-2*x(6)/(91-2*x(6)-1):-x(6) -x(6):-0.1]*x(7);
visualbias=[visualbias visualbias visualbias visualbias];

vb=visualbias(1:15:360);
v=[x(1);x(2)];
Ttheta=[0:pi/12:2*pi-0.001];

Vtheta=Ttheta+vb;
[tx,ty] = pol2cart(Vtheta,1);
T=[tx; ty];
bias_vec=[x(3);x(4)];

distance=sqrt(nansum((T-v).^2)).^x(5);
bias=bias_vec.*distance;
tp=T+bias;
startp=sqrt(nansum(v.^2)).^x(5).*bias_vec;
planp=tp-startp;
[Mtheta,rho] = cart2pol(planp(1,:),planp(2,:));
error_theta=Mtheta-Ttheta;
f=find(error_theta<-pi);
error_theta(f)=error_theta(f)+2*pi;

SE=(error_theta-data2);
SE=SE(:);
Sigm=nanstd(SE);
LLH=log(normpdf( SE,0 ,Sigm));
l2=-nansum(nansum(LLH));


l=l2;




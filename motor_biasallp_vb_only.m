function l=motor_biasallp_vb_only(x,data1)
% TR+TG model

visualbias=[0:x(1) x(1)-0.01:-2*x(1)/(91-2*x(1)-1):-x(1) -x(1):-0.1]*x(2);
visualbias=[visualbias visualbias visualbias visualbias];

vb=visualbias(1:15:360);

error_theta=vb;
f=find(error_theta<-pi);
error_theta(f)=error_theta(f)+2*pi;

SE=(error_theta-data1);
SE=SE(:);
Sigm=nanstd(SE);
LLH=log(normpdf( SE,0 ,Sigm));
l=-nansum(nansum(LLH));







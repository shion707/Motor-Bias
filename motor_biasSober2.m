function l=motor_biasSober2(data1,x)
% Prop Joint


l1=3;
l2=3;
a=x(1);
b=x(2);
biasa=x(3);
biasb=x(4);

x1=cos(a)*l1+cos(b)*l2;
y1=sin(a)*l1+sin(b)*l2;

Ttheta=[0:pi/12:2*pi-0.001];
[tx,ty] = pol2cart(Ttheta,1);

x0=cos(a+biasa)*l1+cos(b+biasb)*l2;
y0=sin(a+biasa)*l1+sin(b+biasb)*l2;
tx=tx+x0;
ty=ty+y0;

for ta=1:24
    dis=sqrt(sum([tx(ta)-x1,ty(ta)-y1].^2));
    x2=x1+(tx(ta)-x1)*1/dis;
    y2=y1+(ty(ta)-y1)*1/dis;
    d2=sqrt(x2.^2+y2.^2);
    cs=(l1^2-l2^2+d2^2)/2/l1/d2;
    theta=acos(cs);
    if x2>0
        deg2=atan(y2/x2);
    else
        deg2=pi-atan(-y2/x2);
    end
    a2=deg2-theta;
    da=a2-a;

    cs3=(l1^2+l2^2-d2^2)/2/l1/l2;
    theta3=acos(cs3);
    b3=pi+a2-theta3;
    db=b3-b;

    rx=cos(a+da+biasa)*l1+cos(b+db+biasb)*l2;
    ry=sin(a+da+biasa)*l1+sin(b+db+biasb)*l2;
    try
        [etheta,~] = cart2pol(rx-x0,ry-y0);
    catch
        etheta=0;
    end
    error_theta(ta)=etheta-Ttheta(ta);
end
f=find(error_theta<-pi);
error_theta(f)=error_theta(f)+2*pi;

SE=(error_theta-data1);
SE=SE(:);
Sigm=nanstd(SE);
LLH=log(normpdf( SE,0 ,Sigm));
l=-nansum(nansum(LLH));


end

 function f = lsc(x,p,t,num_e)
 x=x';
 total=0;
 I=[1,0];
 for i=1:num_e
     Temp=0.5*(p(i)-t(i,1:2)*x)^2+0.5*1*(I*x-1)^2;
     total = total + Temp;
 end
  f=total  ;
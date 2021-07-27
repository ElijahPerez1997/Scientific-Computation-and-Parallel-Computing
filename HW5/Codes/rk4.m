function [x,t] = rk4(f,a,b,N,x0)
x=zeros(length(x0),N);
t=zeros(1,N);

h=(b-a)/N;
t(1)=a;
x(:,1)=x0;

for i=1:N
    k1=h*f(t(i),x(:,i));
    k2=h*f(t(i)+(h/2),x(:,i)+(k1/2));
    k3=h*f(t(i)+(h/2),x(:,i)+(k2/2));
    k4=h*f(t(i)+h,x(:,i)+k3);
    
    x(:,i+1)=x(:,i)+((k1+(2*k2)+(2*k3)+k4)/6);
    t(i+1)=a+i*h;
end
end


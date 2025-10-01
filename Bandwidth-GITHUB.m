% please input the x, y coordinates of the trace. x and y are two vectors.
N=length(x);   %Number of sampling points
Fs=2;   %Sampling frequency

theta1=0;dx=0;dy=0;dxn=0;dyn=0;Y1=0;P1_1=0;P1_2=0;  %Initialization 
for i=1:N-1   %Unit vector
   dx(i)=x(i+1)-x(i);
   dy(i)=y(i+1)-y(i);
   dxn(i)=dx(i)/sqrt(dx(i)*dx(i)+dy(i)*dy(i));  %x component of the unit vector  ni
   dyn(i)=dy(i)/sqrt(dx(i)*dx(i)+dy(i)*dy(i));  %y component of the unit vector
end
for j=2:N-1  %Calculate angle between ni and ni-1
    sint=(dxn(j-1)*dyn(j)-dxn(j)*dyn(j-1)); 
    cost=(dxn(j-1)*dxn(j)+dyn(j-1)*dyn(j));
    if(cost<0 && sint>0)
        theta1(j)=acos(cost); 
    else if(cost<0 && sint<0)
        theta1(j)=-acos(cost);
        else theta1(j)=asin(sint);
        end
    end
  %   theta1(j)= abs(theta1(j)); %You may try what happens using magnitude
end
M=N-1;
m=1:M;
Y1=fft(theta1);
p1_2=abs(Y1/M); 
p1_1=p1_2(1:M/2+1); 
p1_1(2:end-1)=2*p1_1(2:end-1);
f=Fs*(0:(M/2))/M;   %frequencies
f1=f;p1=p1_1; 
plot(f,p1_1,'r-.');  
xlabel('f*'); ylabel('A'); %frequency spectrum
Size_f= size(f); Totalfslot=Size_f(2);

f0=0;
for l=1:Totalfslot  % find the f with maximum A (central frequency)
    if (p1(l)==max(p1))
        f0=f(l);
    end
end

Nzeropoint=0; zeropoints=0;
pbw=p1-max(p1)/2;  %This is called a "-6dB bandwidth"
for k=2:Totalfslot     % find zeropoints
    if pbw(k)*pbw(k-1)<=0
        Nzeropoint=Nzeropoint+1;
        zeropoints(Nzeropoint)=f(k-1)-pbw(k-1)*(f(k)-f(k-1))/(pbw(k)-pbw(k-1));
    end
end

f_low=min(zeropoints); f_high=max(zeropoints);  
lowerbandpass=f0-f_low; %Corner frequencies and lower bandpass

bandwidth1=0; %Sum of bandpasses
if f(1)<=0
  if Nzeropoint==1
      bandwidth1=max(f)-zeropoints(1);
  elseif mod(Nzeropoint,2)==0   %even
    for I=1:Nzeropoint/2
        bandwidth1=bandwidth1+zeropoints(2*I)-zeropoints(2*I-1);
    end
  else
    for I=1:(Nzeropoint-1)/2  %odd
        bandwidth1=bandwidth1+zeropoints(2*I)-zeropoints(2*I-1);
    end
    bandwidth1=bandwidth1+max(f)-zeropoints(Nzeropoint);
  end
else
   if Nzeropoint==1
      bandwidth1=zeropoints(1);
   elseif Nzeropoint==2
      bandwidth1=zeropoints(1)+max(f)-zeropoints(2);   
   elseif mod(Nzeropoint,2)==0   %even >2
      for I=1:Nzeropoint/2-1
          bandwidth1=bandwidth1+zeropoints(2*I+1)-zeropoints(2*I);
      end
      bandwidth1=bandwidth1+zeropoints(1)+max(f)-zeropoints(2); 
   else
    for I=1:(Nzeropoint-1)/2  %odd >2
        bandwidth1=bandwidth1+zeropoints(2*I+1)-zeropoints(2*I)
    end
    bandwidth1=bandwidth1+zeropoints(1);
   end 
end
bandwidth2=0;  %Bandwidth
if f(1)<=0
   if pbw(Totalfslot)<0
       bandwidth2=zeropoints(Nzeropoint)-zeropoints(1);
   else
       bandwidth2=max(f)-zeropoints(1);
   end
else
    if pbw(Totalfslot)<0
        bandwidth2=zeropoints(Nzeropoint);
    else
        bandwidth2=max(f);
    end
end

  
 

 


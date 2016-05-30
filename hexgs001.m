%

clear all

source "functionsHex.m"

   N=6; Gamma=0.4; Delta=0; h=0.92;    

   sigmaz=[1,0;0,-1];

   r_vec=[0.5,sqrt(3)/2;1,0;0.5,-sqrt(3)/2;-0.5,-sqrt(3)/2;-1,0;-0.5,sqrt(3)/2];

   Nk=10;

 qx=(-2.5*pi:(1/Nk):2.5*pi);
 qy=(-2.5*pi:(1/Nk):2.5*pi);
 %qx=2*pi;
 %qy=0;
   [states_array,energies_array]=DiagpbHmat(N,Gamma,Delta,h);
%   SQxx=zeros(length(qx),length(qy));
%   SQyy=zeros(length(qx),length(qy));
%   SQzz=zeros(length(qx),length(qy));  
%   SQxy=zeros(length(qx),length(qy));
%   SQyx=zeros(length(qx),length(qy)); 

       for delta=1:N
          for deltap=1:N 
              thetad=theta(delta,N);
              thetadp=theta(deltap,N);
              sigmax1=sigmaxp(thetad);
              sigmay1=sigmayp(thetad);
              sigmax2=sigmaxp(thetadp);
              sigmay2=sigmayp(thetadp);
              GSCorrxx(delta,deltap)=gsCorr(delta,deltap,sigmax1,sigmax2,N,states_array(:,1));
              GSCorryy(delta,deltap)=gsCorr(delta,deltap,sigmay1,sigmay2,N,states_array(:,1));
              GSCorrzz(delta,deltap)=gsCorr(delta,deltap,sigmaz,sigmaz,N,states_array(:,1));
              GSCorrxy(delta,deltap)=gsCorr(delta,deltap,sigmax1,sigmay2,N,states_array(:,1));
              GSCorryx(delta,deltap)=gsCorr(delta,deltap,sigmay1,sigmax2,N,states_array(:,1));
         end
      end

   for xx=1:length(qx)
   for yy=1:length(qy)

       Sqxx=0;
       Sqyy=0;
       Sqzz=0;
       Sqxy=0;
       Sqyx=0;
       
       for delta=1:N
          for deltap=1:N 
              vec=[qx(xx),qy(yy)]*(r_vec(delta,:)-r_vec(deltap,:))';
              Sqxx=Sqxx+(GSCorrxx(delta,deltap))*exp(i*vec);
              Sqyy=Sqyy+(GSCorryy(delta,deltap))*exp(i*vec);
              Sqzz=Sqzz+(GSCorrzz(delta,deltap))*exp(i*vec);
              Sqxy=Sqxy+(GSCorrxy(delta,deltap))*exp(i*vec);
              Sqyx=Sqyx+(GSCorryx(delta,deltap))*exp(i*vec);
              qxhat=xhat(qx(xx),qy(yy));
              qyhat=xhat(qy(yy),qx(xx));
              SQxx(yy,xx)=(1-qxhat^2)*Sqxx;
              SQyy(yy,xx)=(1-qyhat^2)*Sqyy;
              SQzz(yy,xx)=Sqzz;
              SQxy(yy,xx)=-(qxhat*qyhat)*Sqxy;
              SQyx(yy,xx)=-(qyhat*qxhat)*Sqyx;
          endfor
       endfor

   endfor
   endfor


       
#   data=real(SQxx)+real(SQyy)+real(SQzz)+real(SQxy)+real(SQyx);

#save aN6SQG06h15.dat data    
#save longSqxx002.dat Sqxx

dataxx=real(SQxx);
datayy=real(SQyy);
datazz=real(SQzz);
dataxy=real(SQxy);
datayx=real(SQyx);

save "SQxxh092.dat" dataxx
save "SQyyh092.dat" datayy
save "SQzzh092.dat" datazz
save "SQxyh092.dat" dataxy
save "SQyxh092.dat" datayx



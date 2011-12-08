function [H tau_s] = H_compute(Ish_x,Ish_y,M,ks,idx,tau_s)
% This function implement for H computing
% H(M.*||Delta_Ish||_2-tau_s)
[Ix,Iy] = size(Ish_x);
Delta_Ish = zeros(Ix,Iy);
ratio= 1.1;
if idx==1
   angle1=[];
   angle2=[];
   angle3=[];
   angle4=[];

   for i= ks(idx)+1:Ix-ks(idx)
      for j= ks(idx)+1:Iy-ks(idx)
         theta = atan(Ish_x(i,j)/Ish_y(i,j));
         Delta_Ish(i,j)=sqrt(Ish_x(i,j)^2+Ish_y(i,j)^2);
         if (-pi/2 <= theta) && (theta <= -pi/4)
              angle1=[angle1;Delta_Ish(i,j)];
           
         elseif (-pi/4 < theta) && (theta <= 0)
              angle2=[angle2;Delta_Ish(i,j)];
       
         elseif (0 < theta) && (theta <= pi/4)
              angle3=[angle3;Delta_Ish(i,j)];
           
         else
              angle4=[angle4;Delta_Ish(i,j)];
           
         end
      end
   end

   count=ceil(20*(2*ks(idx)+1));
   tau=[];
   for i=1:4
      switch i
          case 1
              sam=angle1;
          case 2
              sam=angle2;
          case 3
              sam=angle3;
          case 4
              sam=angle4;
      end
      sam=sort(sam,1,'descend');
      tau=[tau,sam(count)];
   end
   tau_s=min(tau);

else
   for i= ks(idx)+1:Ix-ks(idx)
      for j= ks(idx)+1:Iy-ks(idx)
         Delta_Ish(i,j)=sqrt(Ish_x(i,j)^2+Ish_y(i,j)^2);
      end
   end
   tau_s=tau_s/ratio;
end

H=zeros(Ix,Iy);
for i= ks(idx)+1:Ix-ks(idx)
   for j= ks(idx)+1:Iy-ks(idx)
       if (tau_s <= M(i,j)*Delta_Ish(i,j))
           H(i,j)=1;
       end
   end
end
function [M tau_r] = M_compute(Im_x,Im_y,ks,idx,tau_r) 
% This function implement for M computing
% M = H(rmap-tau_r)
[Ix,Iy] = size(Im_x);
rmap = zeros(Ix,Iy);
ksize = 2*ks(idx)+1;
ratio=1.1;
if idx==1
   angle1=[];
   angle2=[];
   angle3=[];
   angle4=[];
   pi=3.1415926;
   for i= ks(idx)+1:Ix-ks(idx)
      
      sumabs_x = 0;
      sumabs_x_col = zeros(1,ksize);
      sum_x = 0;
      sum_x_col = zeros(1,ksize);
      
      sumabs_y = 0;
      sumabs_y_col = zeros(1,ksize);
      sum_y = 0;
      sum_y_col = zeros(1,ksize);
      
      for j= ks(idx)+1:Iy-ks(idx)
          if ks(idx)+1 == j
              for n = j-ks(idx):j+ks(idx)
                 offset = mod(n-1,ksize)+1;
                 for m = i-ks(idx):i+ks(idx)
                     
                    sumabs_x_col(offset) = sumabs_x_col(offset)+abs(Im_x(m,n));
                    sum_x_col(offset) = sum_x_col(offset)+Im_x(m,n);
                    
                    sumabs_y_col(offset) = sumabs_y_col(offset)+abs(Im_y(m,n));
                    sum_y_col(offset) = sum_y_col(offset)+Im_y(m,n);
                    
                 end
              end
              sumabs_x = sum(sumabs_x_col);
              sum_x = sum(sum_x_col);
              sumabs_y = sum(sumabs_y_col);
              sum_y = sum(sum_y_col);
          else
              
              sumabs_x_col_j = 0;
              sum_x_col_j = 0;
              
              sumabs_y_col_j = 0;
              sum_y_col_j = 0;

              for m=i-ks(idx):i+ks(idx)
                  
                 sumabs_x_col_j = sumabs_x_col_j+abs(Im_x(m,j+ks(idx)));
                 sum_x_col_j = sum_x_col_j+Im_x(m,j+ks(idx));

                 sumabs_y_col_j = sumabs_y_col_j+abs(Im_y(m,j+ks(idx)));
                 sum_y_col_j = sum_y_col_j+Im_y(m,j+ks(idx));
                 
              end
              offset = mod(j+ks(idx)-1,ksize)+1;
              
              sumabs_x = sumabs_x-sumabs_x_col(offset)+sumabs_x_col_j;
              sumabs_x_col(offset) = sumabs_x_col_j;
              sum_x = sum_x-sum_x_col(offset)+sum_x_col_j;
              sum_x_col(offset) = sum_x_col_j;
              
              sumabs_y = sumabs_y-sumabs_y_col(offset)+sumabs_y_col_j;
              sumabs_y_col(offset) = sumabs_y_col_j;
              sum_y = sum_y-sum_y_col(offset)+sum_y_col_j;
              sum_y_col(offset) = sum_y_col_j;
          end
          
          %% Sample Rmap computing
%           sumabs_x=0;
%           sum_x=0;
%           sumabs_y=0;
%           sum_y=0;
%           for m=i-ks(idx):i+ks(idx)
%               for n=j-ks(idx):j+ks(idx)
%                   sumabs_x = sumabs_x+abs(Im_x(m,n));
%                   sum_x = sum_x+Im_x(m,n);
%                   sumabs_y = sumabs_y+abs(Im_y(m,n));
%                   sum_y = sum_y+Im_y(m,n);
%               end
%           end
          %%
          rmap_x=sum_x/(sumabs_x+0.5);
          rmap_y=sum_y/(sumabs_y+0.5);
          theta=atan(rmap_x/rmap_y);
          rmap(i,j)=sqrt(rmap_x^2+rmap_y^2);
       
          if (-pi/2 <= theta) && (theta <= -pi/4)
              angle1=[angle1;rmap(i,j)];
           
          elseif (-pi/4 < theta) && (theta <= 0)
              angle2=[angle2;rmap(i,j)];
       
          elseif (0 < theta) && (theta <= pi/4)
              angle3=[angle3;rmap(i,j)];
           
          else
              angle4=[angle4;rmap(i,j)];
           
          end
      end      
   end

   count=ceil(0.5*sqrt((2*ks(idx)+1)^2*Ix*Iy));
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
   tau_r=min(tau);

else
   for i= ks(idx)+1:Ix-ks(idx)
       
      sumabs_x = 0;
      sumabs_x_col = zeros(1,ksize);
      sum_x = 0;
      sum_x_col = zeros(1,ksize);
      
      sumabs_y = 0;
      sumabs_y_col = zeros(1,ksize);
      sum_y = 0;
      sum_y_col = zeros(1,ksize);
       
      for j= ks(idx)+1:Iy-ks(idx)
           if ks(idx)+1==j
              for n = j-ks(idx):j+ks(idx)
                 offset = mod(n-1,ksize)+1;
                 for m = i-ks(idx):i+ks(idx)
                     
                    sumabs_x_col(offset) = sumabs_x_col(offset)+abs(Im_x(m,n));
                    sum_x_col(offset) = sum_x_col(offset)+Im_x(m,n);
                    
                    sumabs_y_col(offset) = sumabs_y_col(offset)+abs(Im_y(m,n));
                    sum_y_col(offset) = sum_y_col(offset)+Im_y(m,n);
                    
                 end
              end
              sumabs_x = sum(sumabs_x_col);
              sum_x = sum(sum_x_col);
              sumabs_y = sum(sumabs_y_col);
              sum_y = sum(sum_y_col);
          else
              
              sumabs_x_col_j = 0;
              sum_x_col_j = 0;
              
              sumabs_y_col_j = 0;
              sum_y_col_j = 0;
              for m=i-ks(idx):i+ks(idx)
                  
                 sumabs_x_col_j = sumabs_x_col_j+abs(Im_x(m,j+ks(idx)));
                 sum_x_col_j = sum_x_col_j+Im_x(m,j+ks(idx));

                 sumabs_y_col_j = sumabs_y_col_j+abs(Im_y(m,j+ks(idx)));
                 sum_y_col_j = sum_y_col_j+Im_y(m,j+ks(idx));
                 
              end
              offset = mod(j+ks(idx)-1,ksize)+1;
              
              sumabs_x = sumabs_x-sumabs_x_col(offset)+sumabs_x_col_j;
              sumabs_x_col(offset) = sumabs_x_col_j;
              sum_x = sum_x-sum_x_col(offset)+sum_x_col_j;
              sum_x_col(offset) = sum_x_col_j;
              
              sumabs_y = sumabs_y-sumabs_y_col(offset)+sumabs_y_col_j;
              sumabs_y_col(offset) = sumabs_y_col_j;
              sum_y = sum_y-sum_y_col(offset)+sum_y_col_j;
              sum_y_col(offset) = sum_y_col_j;
           end

           %% Sample Rmap computing
%            sumabs_x=0;
%            sum_x=0;
%            sumabs_y=0;
%            sum_y=0;
%            for m=i-ks(idx):i+ks(idx)
%                for n=j-ks(idx):j+ks(idx)
%                    sumabs_x = sumabs_x+abs(Im_x(m,n));
%                    sum_x = sum_x+Im_x(m,n);
%                    sumabs_y = sumabs_y+abs(Im_y(m,n));
%                    sum_y = sum_y+Im_y(m,n);
%                end
%            end

           rmap_x=sum_x/(sumabs_x+0.5);
           rmap_y=sum_y/(sumabs_y+0.5);
           rmap(i,j)=sqrt(rmap_x^2+rmap_y^2);
       end
   end
   tau_r=tau_r/ratio;
end

M=zeros(Ix,Iy);

for i= ks(idx)+1:Ix-ks(idx)
   for j= ks(idx)+1:Iy-ks(idx)
      if (tau_r <= rmap(i,j))
          M(i,j)=1;
      end
   end
end
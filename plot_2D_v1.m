function A=plot_2D_v1(w,k)

% get the mask with length 2*k+1 x 2*k+1
% k could be an integer or not an integer
% w is the area under the curve for each bar
%
% Ref. Antonio Cicone, Haomin Zhou. "Multidimensional Iterative Filtering method
%      for the decomposition of high-dimensional non-stationary signals".
%      Cambridge Core in Numerical Mathematics: Theory, Methods and
%      Applications, Volume 10, Issue 2, Pages 278-298, 2017.
%      doi:10.4208/nmtma.2017.s05
%
%      Stefano Sfarra, Antonio Cicone, Bardia Yousefi, Stefano Perilli,
%      Leonardo Robol, Xavier P.V. Maldague.
%      "Maximizing the detection of thermal imprints in civil engineering
%      composites after a thermal stimulus - The contribution of an
%      innovative mathematical pre-processing tool: the 2D Fast Iterative
%      Filtering algorithm. Philosophy, comparisons, numerical, qualitative
%      and quantitative results". 2021. Submitted
% 
m=length(w);

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        A=zeros(k);
        
        N=m/k;
        w=[w zeros(1,m)];
%         for i=1:k+1
%             s=(i-1)*N+1;
%             t=i*N;
%             
%             %s1=s-floor(s);
%             s2=ceil(s)-s;
%             
%             t1=t-floor(t);
%             %t2=ceil(t)-t;
%             
%             if floor(t)<1
%                 disp('Ops')
%             end
%             A(1,i)=sum(w(m+ceil(s):m+floor(t)))+s2*w(m+ceil(s))+t1*w(m+floor(t));
%         end
        for i=1:k+1
            for j=i:k+1
                if i==1 && j==1
                    
                    t=N/2;
                    
                    t1=t-floor(t);
                    
                    if floor(t)<1
                        disp('Ops')
                    end
                    A(i,j)=2*(sum(w(2:floor(t)))+t1*w(floor(t))+w(1));
                else
                    d=sqrt((i-1)^2+(j-1)^2);
                    s=(d-1)*N+1;
                    t=d*N;
                    
                    %s1=s-floor(s);
                    s2=ceil(s)-s;
                    
                    t1=t-floor(t);
                    %t2=ceil(t)-t;
                    
                    if floor(t)<1
                        disp('Ops')
                    end
                    A(i,j)=sum(w(ceil(s):floor(t)))+s2*w(ceil(s))+t1*w(floor(t));
                end
            end
        end
        A=A'+A-diag(diag(A)); % we complete the fourth quadrant         
        A=[rot90(A(2:end,2:end),2) rot90(A(:,2:end),1); rot90(A(2:end,:),3) A]; % Final filter size 2*k+1 x 2*k+1        
%         surf(A)
%         pause(1)
    else   % if the mask length is not an integer
        disp('Need to write the code!')
        A=[];
        return
%         new_k=floor(k);
%         extra = k-new_k;
%         c=(2*m+1)/(2*new_k+1+2*extra);
%         
%         a=zeros(1,2*new_k+3);
%         
%         t=extra*c+1;
%         t1=t-floor(t);
%         %t2=ceil(t)-t;
%         if k<0
%             disp('Ops')
%             a=[];
%             return
%         end
%         a(1)=sum(w(1:floor(t)))+t1*w(floor(t));
%         
%         for i=2:2*new_k+2
%             s=extra*c+(i-2)*c+1;
%             t=extra*c+(i-1)*c;
%             %s1=s-floor(s);
%             s2=ceil(s)-s;
%             
%             t1=t-floor(t);
%             
%             
%             a(i)=sum(w(ceil(s):floor(t)))+s2*w(ceil(s))+t1*w(floor(t));
%         end
%         t2=ceil(t)-t;
%         
%         a(2*new_k+3)=sum(w(ceil(t):n))+t2*w(ceil(t));
    end
 else % We need a filter with more points than MM, we use interpolation
     disp('Need to write the code!')
        A=[];
        return
%     dx=0.01;
%     % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
%     % filter of length 62*2 in the physical space
%     f=w/dx; % function we need to interpolate
%     dy=m*dx/k;
%     b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
%     if size(b,1)>size(b,2)
%         b=b.';
%     end
%     if size(b,1)>1
%         fprintf('\n\nError!')
%         disp('The provided mask is not a vector!!')
%         a=[];
%         return
%     end
%     a=[fliplr(b(2:end)) b]*dy;
%     if abs(norm(a,1)-1)>10^-14
%         %         fprintf('\n\nError\n\n')
%         %         fprintf('Area under the mask = %2.20f\n',norm(a,1))
%         %         fprintf('it should be equal to 1\nWe rescale it using its norm 1\n\n')
%         a=a/norm(a,1);
%     end
end

end
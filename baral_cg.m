function x = baral_cg(N,n_iter)
%   Conjugate Gradient Method.
% Part of Coursework for Mathmatical modeling at the
% University of Texas at El Paso
% Submitted by: Pratik Baral on May 6/2017
% Submitted to: Dr. Natasha Sharma
% How to Run

%Example run as conjgrad(5, 50) 
% 5 is the dimension of the square matrix A and 50 is the no of iteration
% Note the input here is N = dimension of the square matrix and default 
%tolerance is 1e-7 and h = .001

 tol=1e-7;
A = zeros(N,N);
b = ones(N,1);
h = .001;
for i = 1:N
    
   for j = 1: N
%        a[i][j]=0;
       if i==j
            A(i,j)=4/h^2;
       end
       if int8(i/j)==int8(j/i) & abs(i-j)==1
           A(i,j)=-1/h^2;
       end
       if abs(i-j)==N
           A(i,j)=-1/h^2;
       end           
   end    
end
  
       
  
    x = zeros(N,1);
    r = b - A*x;
    if norm(r) < tol
        return
    end
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
   fprintf('i \t norm \n ')
    for k = 1:n_iter;
       r = r - t*z;
      
       fprintf('%d \t %f\n ', k, norm(r))
       if( norm(r) < tol )
            return;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
    end
end
 
% reference https://www.mathworks.com/matlabcentral/fileexchange/22494-conjugate-gradient-method
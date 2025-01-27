%====================================================================
% this function is the same as fvnew1_revised_agian.m except that 
% inv(A) is computed by interaitve methods, and with recycling technique
% whereas, in fvnew1_revised_agian.m inv(A) is computed by LU
% decomposition.
% the parametric system there is(a second order system):
% (M1+dM2)d^2x/dt^2+theta(D1+dD2)dx/dt+(T1+1/dT2+dT3)x=B
%                                 y(t)  =Cx
% where matrices M1=Bmass1,M2=Bmass2,D1=Bdamp1,D2=Bdamp2,T1=Bstiff1,T2=Bstiff2,T3=Bstiff3,B=Bload are stored in orig_matrices.mat
%======================================================================

% Kapil
% BiCG and BiCGSTAB Code
addpath('../');

n=size(B,1);

p=4;   %number of parameters


Et=sparse(n,p*n);

Et(:,1:n)=E0;

Et(:,n+1:2*n)=E1;

Et(:,2*n+1:3*n)=A1;

Et(:,3*n+1:4*n)=A2;


%===================================================
% compute the number  of columns s included in V
%===================================================


V=sparse(n,1);   %initialize V

s1=0;  %the number of columns in V computed from the 0th to j-1th order moments

s=0;   %the number of columns in V computed from the oth to the current jth order moments


%-----------------------------------------------------------
%   compute first column in V by the 0th moment
%-----------------------------------------------------------

% Kapil  
m = 25;
k = 20;
tol = 1e-8;
maxit  = 100;
c = ones(n,1);
[Lp,Up] = luinc(coeff, 1e-4);
c1 = Up'\c;

num_systems = num_systems + 1;

% Direct
b=coeff\B;
 
% No recycling
tic
B1 = Lp\B;         % Central Preconditioning; 
                   % RBiCSTAB does not apply Lp
[b1, ~, iter, ~] = rbicgstab2(coeff, B1, c1, tol, maxit, x01, x0_tilde1, ...
    [], [], [], [], [], [], Lp, Up);
x01 = b1;
b1 = Up\b1;        % Central Preconditioning; 
                   % RBiCSTAB does not apply Up and return
vec_bicgstab_time(num_systems)= toc;
vec_bicgstab_iter(num_systems) = iter;
vec_bicgstab_res(num_systems) = norm(B - coeff*b1)/norm(B);

% With recycling
tic
B2 = Lp\B;          % Central Preconditioning; 
                    % RBiCSTAB does not apply Lp
[b2, ~, ~, iter, ~, ~, U1, U1_tilde, C1, C1_tilde, C1_hat, C1_tilde_hat] = ...
    rbicg2(coeff, B2, c1, tol, maxit, x02, x0_tilde2, U1, U1_tilde, 1, m, k, Lp, Up);
x02 = b2;
b2 = Up\b2;         % Central Preconditioning; 
                    % RBiCSTAB does not apply Up and return
vec_recycling_time(num_systems)= toc;
vec_recycling_iter(num_systems) = iter;
vec_recycling_res(num_systems) = norm(B - coeff*b2)/norm(B);
    
%-------------------------------------------------------------------

 V(:,1)=b/norm(b);

 s=s+1;

%==================================================================
% compute the following columns in V starting from the 1th moment
%==================================================================

for i=1:JJ    % ith moment, maximum to the JJth order moments
    
    i
        
    s2=s;      % remember the number of columns in V corresponding to 
               % the 0th to the i-1th order moments
    
       
   %--------------------------------------------------------------------------------------
   % compute the p blocks  in the ith order moment, each block
   % corresponding to the coefficient matrix M_t before the each parameter s_t
   % t=1,2,...,p
   %--------------------------------------------------------------------------------------
   
   s0=s;
   for t=1:p
       
       t
       
       if s1==s0            % this means all columns which are computed from
                            % the i-1th order moments are linearly dependent with
          break             % precious columns, they are deleted imediately when they are
                            % detected as almost zeros, there is no columns
       else                 % corresponding to the i-1th order moments, we loose the
                            % basis to continue the algorithm to compute the ith
                            % order moments, which need i-1th order moments to
                             % multiply the coefficient matrix, therefore,the algorithm stops.
          for j=s1+1:s0   
             
              if t==1
                  
                   w0=1e9*Et(:,(t-1)*n+1:t*n)*V(:,j);      % w0=M_t*V(:,j); 
                  
              elseif t==3
                  
                   w0=1e4*Et(:,(t-1)*n+1:t*n)*V(:,j);      % w0=M_t*V(:,j); 
                   
              else
                  
                    w0=1e15*Et(:,(t-1)*n+1:t*n)*V(:,j); 
                    
              end

              if norm(w0)>1e-6
                
                %==========================================================
                %---------------------------------------------------------
                
                % Kapil
                num_systems = num_systems + 1;
                
                 % Direct
                 w=coeff\w0;
                
%                  % Iterative with BiCGSTAB
%                  tic
%                  [w1,~,~,~,~] = bicgstab(coeff,w0,tol,maxit,Lp,Up);
%                  Initial Guess???
%                  toc
%                  norm(w0 - coeff*w1)/norm(w0) 
%                  % keyboard
                 
                 % Iterative with BiCGSTAB BUT my code because Matlab
                 % does right preconditioning. We do mixed precond.
                 tic
                 w01 = Lp\w0;       % Central Preconditioning; 
                                    % RBiCSTAB does not apply Lp
                 [w1, ~, iter, ~] = rbicgstab2(coeff, w01, c1, tol, maxit, ...
                     x01, x0_tilde1, [], [], [], [], [], [], Lp, Up);
                 x01 = w1;
                 w1 = Up\w1;        % Central Preconditioning; 
                                    % RBiCSTAB does not apply Up and return
                 vec_bicgstab_time(num_systems)= toc;
                 vec_bicgstab_iter(num_systems) = iter;
                 vec_bicgstab_res(num_systems) = norm(w0 - coeff*w1)/norm(w0);
                 
                % Iterative with RBiCG and RBiCGSTAB
                tic
                w02 = Lp\w0;        % Central Preconditioning; 
                                    % RBiCSTAB does not apply Lp
                % Not updating x0_tilde2 so that consistent with w/o recyc
                [w2, ~, iter, ~] = rbicgstab2(coeff, w02, c1, tol, maxit, ...
                    x02, x0_tilde2, U1, U1_tilde, C1, C1_tilde, C1_hat, C1_tilde_hat, Lp, Up);
                x02 = w2;
                w2 = Up\w2;         % Central Preconditioning; 
                                    % RBiCSTAB does not apply Up and return
                vec_recycling_time(num_systems)= toc;
                vec_recycling_iter(num_systems) = iter;
                vec_recycling_res(num_systems) = norm(w0 - coeff*w2)/norm(w0);
                 
                 % flag
                 % relres
                 % iter
                 % plot(resvec)
                 % norm(w-w1)/norm(w)
 
                %=========================================================
                
                 col=s+1;     %the current column in V
     
                 for kk=1:col-1     % orthogonalize the current column in
                                    % V with previous columns
         
                     h=V(:,kk)'*w;
         
                     w=w-h*V(:,kk);
         
                 end
      
                     ww=norm(w);
        
                 if ww>1e-5
            
                    V(:,col)=w/norm(w);
           
                    s=s+1          
           
                 end
                  
              end
              
              
          end
          
       end
   
   end
   
  s1=s2;
  
end
   

     
    

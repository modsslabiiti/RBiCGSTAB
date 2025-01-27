% this function orthognalize the colums of V and delet those linearly depedent colums, Vend is the resulting 

% linearly independent colums and orthognal.
% orthogonalize the colums of V


function Vend=forthognalize(V,tol)

[n1,n2]=size(V);

V0=sparse(n1,n2);

w=V(:,1);

nw=norm(w);

V0(:,1)=w/nw;

s=1;
for j=2:n2;
    
       
    w=V(:,j);
    
    for i=1:s
        
        rr=V0(:,i)'*w;
        
        w=w-rr*V0(:,i);
        
    end
    
    nw=norm(w);
    
    if nw>tol
        
        s=s+1;
        
       V0(:,s)=w/norm(w);
  
    end
  
   
end

Vend=V0(:,1:s);





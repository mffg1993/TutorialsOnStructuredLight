function H = NHermite(n,X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hermite Polynomials 
%   Returns the evaluation on the domain X of the n-th order Hermite
%   Polynomial, where n is a nonnegative interger
%   
%    Manuel Ferrer, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Hn1=ones(size(X));
H=2.*X; 

if n<0  
    error('The index must be 0 or positive')
elseif n==0
    H=Hn1;
elseif n==1
    H=H;
else 
    for nn=2:n
        Hn=2.*X.*H-2.*(nn-1).*Hn1; 
        Hn1=H;
        H=Hn;       
    end
end
end


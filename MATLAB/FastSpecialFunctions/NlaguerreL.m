function LL = NlaguerreL(n,a,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Associated Laguerre Polynomials 
%   Returns the evaluation on the domain X of the associate Laguerre 
%   Polynomial of order n, where n is a nonnegative interger
%   
%    Manuel Ferrer, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL=0;
for m=0:n
LL=LL+((-1)^m).*(factorial(n+a))/(factorial(n-m)*factorial(a+m)*factorial(m)).*(X.^m);
end
end
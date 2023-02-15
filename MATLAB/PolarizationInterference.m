clear all 
close all 
clc 

% Half Simulation window length
xmax=2; 

% Simulation space
x= linspace(-xmax,xmax,100); 
[xx,yy] = meshgrid(x,x); 

%
% Function of a LGB at z=0 and radial order=0 
LGB= @(X,Y,w0,L) sqrt(2/(pi*factorial(abs(L))))*(1/w0)*exp (-(X.^2+Y.^2)/w0^2).* ((sqrt(2)*sqrt (X.^2+Y.^2)./w0).^(abs(L))).*exp(1i*(L*atan2(Y,X)));


% I'm representing a polarized field using a N x N x 2 array. I worked
% everything in the XY (x/horizontal and y/vertical). The arry [:,:,1] is
% the x-polarized field. 

% Here I defined a LGB beam with topological charge L and the polarization
% can be controlled using alpha and beta for controlling the relative
% intensity and phase, respectively 
Field= @(alpha,beta,L) cat(3,cos (alpha).*LGB(xx,yy,1,L),sin(alpha)*exp (1i*beta).*LGB(xx,yy,1,L)); 


%% Example 1 - x polarized beam LGB l=1

EH=Field(0,0,1); 

% Stokes parameters
SH=StokesParameters(EH); 


figure(1)
ellipsemap(SH,20)
title('X-Polarized')


%% Example 2 - R polarized beam LGB l=1

ER=Field(pi/4,pi/2,1); 

% Stokes parameters
SR=StokesParameters(ER); 


figure(2)
title('R-Polarized')
ellipsemap(SR,20)

%% Example 3 - Full Poincare beam

EFP=Field(pi/4,pi/2,0)+Field(pi/4,-pi/2,1); 

% Stokes parameters
SFP=StokesParameters(EFP); 


figure(3)
title('Star Full-Poincare')
ellipsemap(SFP,20)


%% Functions

% Intensity calculation considering a polarized beam
function A=Intensity(EE)
A=abs (EE(:,:,1)).^2+abs (EE(:,:,2)).^2; 
end


% Stokes Parameters considering a polarized beam
function S=StokesParameters(EE)
S0=abs (EE(:,:,1)).^2+abs (EE(:,:,2)).^2; 
S1=abs (EE(:,:,1)).^2-abs (EE(:,:,2)).^2; 
S2=2*real(conj (EE(:,:,2)).*EE(:,:,1)); 
S3=-2*imag(conj (EE(:,:,2)).*EE(:,:,1)); 

S=cat(3,S0,S1,S2,S3); 
end


function ellipsemap(StokesVector,NE)
% New version of Polarization Ellipse map 
% Ferrer, Lopex 2023
% 
% Shows and calculates the polarization ellipse of a field when in
% differentes points of a rectangular grid. The ellipses are calculated
% from the stoke parameters measurements.
%
%
% Parameters: 
%
% S   - Vector that contains the stokes parameters measurements. All of them should
%     have the same dimensions
% NE  - Number of ellipses  NE>=2 
% 
% Other parameters
% 
% fill     - The size of ellipse within its unit cell. f=1 is the max
%               vale
% tresh    - Intesitiy treshhold. If the intensity at r(x,y) is less that
%            S0*tresh, the ellipse is not generated.


%Fill Factor
fill = 0.7;

% Threshhold
thresh=0.18;
% Array size a-> rows, b-> columns
[a,b] = size(StokesVector(:,:,1)); 

S0=StokesVector(:,:,1);

% Stokes parameter normalization
s1=StokesVector (:,:,2)./StokesVector(:,:,1);
s2=StokesVector (:,:,3)./StokesVector(:,:,1);
s3=StokesVector (:,:,4)./StokesVector(:,:,1);
% Angular coordinate - parametrization
t=linspace(0,2*pi,100);
% Points in the grid
puntos_a=round(linspace(1,a,NE*2+1));
puntos_b=round(linspace(1,b,NE*2+1));
% Scaler
Emax=min([round(fill*a/(NE*2+1)),round(fill*b/(NE*2+1))]);
% Intesity plot
imagesc(StokesVector(:,:,1)); colormap gray; set(gca,'YDir','normal');
hold on
% Ellipse generation
for col=puntos_b(2:2:end)
    for ren=puntos_a(2:2:end)
        if S0(ren,col)<max(max(S0))*thresh
            xp=0; yp=0;
        else
            beta=(0.5)*atan2(s3(ren,col),sqrt(s1 (ren,col)^2+s2 (ren,col)^2));
            Emin=tan(beta)*Emax;
            alpha=(0.5)*atan2(s2(ren,col),s1(ren,col));
            xp=Emax*cos(alpha)*cos(t)-Emin*sin(alpha)*sin(t)+col;
            yp=Emax*sin(alpha)*cos(t)+Emin*cos(alpha)*sin(t)-ren+a+1;
            if abs(s3(ren,col)) <10*eps
                plot(xp,yp,'b','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880])
            else
                if s3(ren,col)> 0
                    plot(xp,yp,'b','LineWidth',1.5,'Color','r') % right circular
                else
                    plot(xp,yp,'b','LineWidth',1.5,'Color','b') % left circular
                end
            end
        end
    end
end
hold off; axis equal; axis off;

end 

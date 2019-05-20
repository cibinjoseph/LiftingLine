% Lifting Line theory for fixed flat plate wing
clc; clear;

% Parameters
span = 8;
chord = 1;
ns = 24;
theta = 5*(pi/180);
cA = [100 ; 0 ; 0];
cD = [100 ; 0 ; 0];
vinf = [10 ; 0 ; 0];
rho = 1.2;

% Create geometry
yvec = linspace(0,span,ns+1);
dy = yvec(2)-yvec(1);

cB = zeros(3,ns); cB(1,:) = 0.25*chord;
cC = cB;
cB(2,:) = yvec(1:end-1);
cC(2,:) = yvec(2:end);

cp = zeros(3,ns);
cp(1,:) = 0.75*chord;
cp(2,:) = (cB(2,:)+cC(2,:))*0.5;

ncap = [sin(theta) ; 0; cos(theta)];

% % Plot wing and cp
%plot(cB(1,:),cB(2,:),'ro')
%hold on;
%plot(cC(1,:),cC(2,:),'b*')
%plot(cp(1,:),cp(2,:),'bo')

% AIC
AIC = zeros(ns,ns);
for is = 1:ns   % panel
  cA(2) = cB(2,is);
  cD(2) = cC(2,is);
  for icp = 1:ns   % cp
    AIC(is,icp) = dot(hshoe(cp(:,icp),cA,cB(:,is),cC(:,is),cD,1),ncap);
  end
end

% RHS
RHS = ones(ns,1)*dot(-vinf,ncap);

% Gam
gam = inv(AIC)*RHS;

% Lift
lift = rho*norm(vinf)*gam*dy;

% CL
cl_liftingline = sum(lift)/(0.5*rho*(norm(vinf))^2*span*chord)
cl_2d = 2*pi*theta

function [ attackAngle, stall, colName ] = isstall( controls, xstruct )

%   Calculate angle of attack for each of LRAUV fins and determain if the
%   vehicle is stalling.
%
%   http://en.wikipedia.org/wiki/Angle_of_attack    
%   Based on code written by Rob McEwen: robsFins.m
%   Last modified Jan 07, 2015
%   Ben Raanan


% Initialize global variables
%---------------------------------------------------------------------
% vehicle_coeffs ;

% global Sfin ARe dCL CDc

% rho     =   1025;             % kg/m3
% ARe     =   6.500000;         % n/a     Fin aspect ratio
dCL     =   4.130000;         % n/a     Coef. of Lift Slope
CDc     =   0.600000;         % n/a     Crossflow Coef. of Drag !!try 0.6!!
Cd0     =   0.030000;         % n/a     Min reaction drag coeff
ec      =   0.9;
Sfin    =   1.15e-2;          % m^2   Fin area
bfin    =   18.57e-2;         % m     Fin span
ARe     = 2*((bfin^2)/Sfin);  % n/a   Effective aspect ratio
cr      =  0.0854;            % m     Root chord
ct      =  0.0854;            % m     Tip chord
cm      =  (cr+ct)/2;         % m     Mean chord   


stall_angle    =   15.0*pi/180; % rad

% Position of fins from vehicle center:
colName = {'lowerRud','portElev','upperRud','stbdElev'};
% rowName = {'X','Y','Z'};
xi = 1 ; yi = 2 ; zi = 3 ; % dimention indecies
finPos  = [ -0.633, -0.633, -0.633, -0.633 ;
             0.012, -0.152,  0.012,  0.152 ;
            -0.152,  0.000,  0.152,  0.000 ] ;  % m

% Get and check state variables and control inputs
%---------------------------------------------------------------------
% Get state variables
u  = xstruct.u ; v  = xstruct.v ; w  = xstruct.w ;
p  = xstruct.p ; q  = xstruct.q ; r  = xstruct.r;

if isdeg(xstruct.theta)
    p = deg2rad(p); q = deg2rad(q); r = deg2rad(r);
end

% Get control inputs
if isdeg(controls(:,1))
    elev_ang = deg2rad(controls(:,1)) ; rud_ang = deg2rad(controls(:,2)) ;
else
    elev_ang = controls(:,1) ; rud_ang = controls(:,2) ;
end

% Get temp and sal
T = xstruct.Temp ; uT = 'K' ;
S = xstruct.Sal ;  uS = 'ppt';
rho = SW_Density(T,uT,S,uS);
nu = SW_Kviscosity(T,uT,S,uS);


% Initialize elements of coordinate system transform matrix
%---------------------------------------------------------------------
% s1 = sin(rud_ang); s2 = sin(elev_ang);
% c1 = cos(rud_ang); c2 = cos(elev_ang);


% Calculate fin velocities in body coords
%---------------------------------------------------------------------
v1  = [ u + q*finPos(zi,1) - r*finPos(yi,1) ; % 'lowerRud'
        v - p*finPos(zi,1) + r*finPos(xi,1) ;
        w + p*finPos(yi,1) - q*finPos(xi,1) ] ;

v2  = [ u + q*finPos(zi,2) - r*finPos(yi,2) ; % 'portElev'
        v - p*finPos(zi,2) + r*finPos(xi,2) ;
        w + p*finPos(yi,2) - q*finPos(xi,2) ] ;

v3 	= [ u + q*finPos(zi,3) - r*finPos(yi,3) ; % 'upperRud'
        v + r*finPos(xi,3) - p*finPos(zi,3) ;
        w + p*finPos(yi,3) - q*finPos(xi,3) ] ;

v4  = [ u + q*finPos(zi,4) - r*finPos(yi,4) ; % 'stbdElev'
        v - p*finPos(zi,4) + r*finPos(xi,4) ;
        w + p*finPos(yi,4) - q*finPos(xi,4) ] ;


% Now get angle of attack for each fin
%---------------------------------------------------------------------
norm_v1 = sqrt( v1(1,:).*v1(1,:) + v1(2,:).*v1(2,:) + v1(3,:).*v1(3,:) );
norm_v2 = sqrt( v2(1,:).*v2(1,:) + v2(2,:).*v2(2,:) + v2(3,:).*v2(3,:) );
norm_v3 = sqrt( v3(1,:).*v3(1,:) + v3(2,:).*v3(2,:) + v3(3,:).*v3(3,:) );
norm_v4 = sqrt( v4(1,:).*v4(1,:) + v4(2,:).*v4(2,:) + v4(3,:).*v4(3,:) );

% Reynolds number time-series
norm_v  = [norm_v1; norm_v2; norm_v3; norm_v4];
Re = NaN(size(norm_v));
for k=1:size(norm_v,1)
    Re(k,:) = (norm_v(k,:).*cm)./nu;
end


alpha1 =  rud_ang'  - v1(2,:)./v1(1,:);     % 'lowerRud'
alpha2 =  elev_ang' + v2(3,:)./v2(1,:);     % 'portElev'
alpha3 =  rud_ang'  - v3(2,:)./v3(1,:);     % 'upperRud'
alpha4 =  elev_ang' + v4(3,:)./v4(1,:);     % 'stbdElev'

alpha1( norm_v1 < 0.001 ) = 0.0; 
alpha2( norm_v2 < 0.001 ) = 0.0; 
alpha3( norm_v3 < 0.001 ) = 0.0; 
alpha4( norm_v4 < 0.001 ) = 0.0; 

% Lift coefficients in construction...
%{
% lift and coefficients */
CDC = CDc/ARe;

CL1 = dCL.*alpha1 + CDC*alpha1.*abs(alpha1);
CL2 = dCL.*alpha2 + CDC*alpha2.*abs(alpha2);
CL3 = dCL.*alpha3 + CDC*alpha3.*abs(alpha3);
CL4 = dCL.*alpha4 + CDC*alpha4.*abs(alpha4);

% Note that if stall angle is exceeded: NO LIFT */
% CL1(abs(alpha1) < stall_angle) = 0.0;
% CL2(abs(alpha2) < stall_angle) = 0.0;
% CL3(abs(alpha3) < stall_angle) = 0.0;
% CL4(abs(alpha4) < stall_angle) = 0.0;

aa = 1.0/(pi*ARe*ec);

CD1 = Cd0 + aa.*CL1.*CL1;
CD2 = Cd0 + aa.*CL2.*CL2;
CD3 = Cd0 + aa.*CL3.*CL3;
CD4 = Cd0 + aa.*CL4.*CL4;


% lift and drag forces, in flow coords... */
%---------------------------------------------------------------------
cons = (rho*Sfin)/2.0;

LW1 = cons.*norm_v1.*norm_v1.*CL1;      % positive when the lift vector is close to normal vector */
LW2 = cons.*norm_v2.*norm_v2.*CL2;
LW3 = cons.*norm_v3.*norm_v3.*CL3;
LW4 = cons.*norm_v4.*norm_v4.*CL4;


DW1 = cons.*norm_v1.*norm_v1.*CD1;      % always positive */
DW2 = cons.*norm_v2.*norm_v2.*CD2;
DW3 = cons.*norm_v3.*norm_v3.*CD3;
DW4 = cons.*norm_v4.*norm_v4.*CD4;


LF1 = LW1.*cos(alpha1) + DW1.*sin(alpha1);        % force in the fin normal direction */
LF2 = LW2.*cos(alpha2) + DW2.*sin(alpha2);
LF3 = LW3.*cos(alpha3) + DW3.*sin(alpha3);
LF4 = LW4.*cos(alpha4) + DW4.*sin(alpha4);




figure;
plot(time, LF2)
dynamicDateTicks(gca)
grid on
%}

% rearrange for output
%---------------------------------------------------------------------
attackAngle(:,1) = alpha1 ; % 'lowerRud'
attackAngle(:,2) = alpha2 ; % 'portElev'
attackAngle(:,3) = alpha3 ; % 'upperRud'
attackAngle(:,4) = alpha4 ; % 'stbdElev'


stall(:,1) = abs(alpha1) > stall_angle ; 
stall(:,2) = abs(alpha2) > stall_angle ;
stall(:,3) = abs(alpha3) > stall_angle ;
stall(:,4) = abs(alpha4) > stall_angle ;
stall = double(stall); % logical to mat


% NaN when vehicle is at the surface or when we don't have controller data
%---------------------------------------------------------------------
ind1 = ( xstruct.z<1 | isnan(xstruct.z) )'; ind2=isnan(elev_ang);
ind = ind1==1 | ind2==1;
for k=1:4;
    stall(ind,k) = NaN;
    attackAngle(ind,k) = NaN;
end

end


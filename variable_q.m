clear all; clc; close all;
Lx = 5.0; Ly = 5.0;                                                  % domain size 
gx = 0.0; gy = -0.0; rho1 =2; rho2 =.1; mu1 = 0.5;mu2 = 0.1; gamma = 20;    % parameters
unorth = 0; usouth = 0; veast = 0; vwest = 0;                        % boundary conditions
rad = 0.5; xc = 2.5; yc = 2.5;                                      % initial drop size and location
D = 0.01;
%k1=.001; 
k2=.001; %conductivities in gas and liquid
kg=k2;
%kl=k1;
hlg=10.0; Cpg=1.0;Cpl=1.0;
rcp1=rho1*Cpl;rcp2=rho2*Cpg;
Tsat=1;Tw=10;Tin = 0.5;
Twvec=[5];
%Uncomment the following for Tw variation
Twvec=[10,5,2];
klvec=[.001];
%Uncomment the following for kl variation
%klvec=[0.001,.0025,0.005];

for iT=1:length(Twvec)
    for ik=1:length(klvec)
    Tw=Twvec(iT);
    kl=klvec(ik);
    k1=kl;
    
time = 0.0; plot_freq = 10;

nx = 128; ny = 128; dx = Lx/nx; dy = Ly/ny; dt = 0.0001;
nstep = 300; maxit = 200; maxError = 0.001; omg = 1.5; Nf = 100;

u=zeros(nx+1,ny+2); ut = u   ; uplot = zeros(nx+1,ny+1);
v=zeros(nx+2,ny+1); vt = u   ; vplot = zeros(nx+1,ny+1);
p=zeros(nx+2,ny+2); tmp1 = p ; tmp2  = p; r = p; chi = p; cond=p;rcp=p;
c=zeros(nx+2,ny+2)+Tw; qdot=zeros(nx+2,ny+2); cold=c;  

xf=zeros(1,Nf+2); yf=zeros(1,Nf+2); 
uf=zeros(1,Nf+2); vf=zeros(1,Nf+2); 
qf=zeros(1,Nf+2);

xh = linspace(0,Lx,nx+1)         ; yh = linspace(0,Ly,ny+1);                % velocity points
x  = linspace(-dx/2,Lx+dx/2,nx+2); y  = linspace(-dy/2,Ly+dy/2,ny+2);       % pressure points

r = zeros(nx+2,ny+2) + rho2;                                                % initial density 
mu = zeros(nx+2,ny+2) + mu2;
for i = 2:nx+1; for j = 2:ny+1                   
  if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2)
      r(i,j) = rho1; 
      mu(i,j) = mu1;
      chi(i,j)=1.0;
  end; 
  if((x(i)-xc)^2+(y(j)-yc)^2 < (rad)^2)
      c(i,j)=Tin;
  end;
end; end;
                                          
for l=1:Nf+2; xf(l) = xc - rad*sin(2*pi*(l-1)/Nf); yf(l) = yc + rad*cos(2*pi*(l-1)/Nf); end   % initialize front                                  
            
for is=1:nstep
    d_n(is,iT,ik)=sqrt(4/pi*sum(sum(chi))*Lx*Ly/(nx*ny));
   

% find q_dot on the interface
 for l=2:Nf+1
     
     n1 = -(yf(l+1)-yf(l-1))/sqrt((xf(l+1)-xf(l-1))^2+(yf(l+1)-yf(l-1))^2);
     n2 = (xf(l+1)-xf(l-1))/sqrt((xf(l+1)-xf(l-1))^2+(yf(l+1)-yf(l-1))^2);
     %point is at a distance 1.5dx from x1,y1. 1.5 is an arbitrary param between 1 and 2
     %point outside the drop
     xplus = xf(l)+n1*1.5*dx; 
     yplus = yf(l)+n2*1.5*dx; 
     %point inside the drop
     xmin = xf(l)-n1*1.5*dx; 
     ymin = yf(l)-n2*1.5*dx;
     
     %Calculate Tg, ufg, vfg by bilinear interpolation
     ip = floor((xplus+0.5*dx)/dx)+1; jp = floor((yplus+0.5*dy)/dy)+1;
     ax = (xplus+0.5*dx)/dx-ip+1; ay = (yplus+0.5*dy)/dy-jp+1;	   
     Tg = (1.0-ax)*(1.0-ay)*c(ip,jp) + ax*(1.0-ay)*c(ip+1,jp) + (1.0-ax)*ay*c(ip,jp+1) + ax*ay*c(ip+1,jp+1);
     ip = floor(xplus/dx)+1; jp = floor((yplus+0.5*dy)/dy)+1;
     ax = xplus/dx-ip+1; ay = (yplus+0.5*dy)/dy-jp+1;	   
     ufg = (1.0-ax)*(1.0-ay)*u(ip,jp) + ax*(1.0-ay)*u(ip+1,jp) + (1.0-ax)*ay*u(ip,jp+1) + ax*ay*u(ip+1,jp+1);
     ip = floor((xplus+0.5*dx)/dx)+1; jp = floor(yplus/dy)+1;
     ax = (xplus+0.5*dx)/dx-ip+1; ay = yplus/dy-jp+1;
     vfg = (1.0-ax)*(1.0-ay)*v(ip,jp) + ax*(1.0-ay)*v(ip+1,jp) + (1.0-ax)*ay*v(ip,jp+1) + ax*ay*v(ip+1,jp+1);
   
     %Calculate Tl, ufl, vfl by bilinear interpolation
     ip = floor((xmin+0.5*dx)/dx)+1; jp = floor((ymin+0.5*dy)/dy)+1;
     ax = (xmin+0.5*dx)/dx-ip+1; ay = (ymin+0.5*dy)/dy-jp+1;	   
     Tl = (1.0-ax)*(1.0-ay)*c(ip,jp) + ax*(1.0-ay)*c(ip+1,jp) + (1.0-ax)*ay*c(ip,jp+1) + ax*ay*c(ip+1,jp+1);
     ip = floor(xmin/dx)+1; jp = floor((ymin+0.5*dy)/dy)+1;
     ax = xmin/dx-ip+1; ay = (ymin+0.5*dy)/dy-jp+1;	   
     ufl = (1.0-ax)*(1.0-ay)*u(ip,jp) + ax*(1.0-ay)*u(ip+1,jp) + (1.0-ax)*ay*u(ip,jp+1) + ax*ay*u(ip+1,jp+1);
     ip = floor((xmin+0.5*dx)/dx)+1; jp = floor(ymin/dy)+1;
     ax = (xmin+0.5*dx)/dx-ip+1; ay = ymin/dy-jp+1;
     vfl = (1.0-ax)*(1.0-ay)*v(ip,jp) + ax*(1.0-ay)*v(ip+1,jp) + (1.0-ax)*ay*v(ip,jp+1) + ax*ay*v(ip+1,jp+1);
   
     %Calculate /dot{q} on the interface assuming Tsat at interface
     qf(l) = (kg*(Tg-Tsat) - kl*(Tsat-Tl))/(1.5*dx); %distance 1.5dx from x1,y1
     %qf(l)=1.0;
     Vn = (ufg+ufl)*n1/2.0 + (vfg+vfl)*n2/2.0 - 0.5*qf(l)/hlg*(1/rho1+1/rho2);
     uf(l)=Vn*n1;
     vf(l)=Vn*n2;
     %Smooth
     ds=sqrt((xf(l+1)-xf(l))^2+(yf(l+1)-yf(l))^2);
     qf1 = qf(l)*ds;
     ip = floor((xf(l)+0.5*dx)/dx)+1; jp = floor((yf(l)+0.5*dy)/dy)+1;
     ax = (xf(l)+0.5*dx)/dx-ip+1;     ay = (yf(l)+0.5*dy)/dy-jp+1;	  
     qdot(ip,jp)     = qdot(ip,jp)     + (1.0-ax)*(1.0-ay)*qf1/dx/dy;
     qdot(ip+1,jp)   = qdot(ip+1,jp)   + ax*(1.0-ay)*qf1/dx/dy;
     qdot(ip,jp+1)   = qdot(ip,jp+1)   + (1.0-ax)*ay*qf1/dx/dy;
     qdot(ip+1,jp+1) = qdot(ip+1,jp+1) + ax*ay*qf1/dx/dy; 
 end;

 xf(2:Nf+1) = xf(2:Nf+1) + dt*uf(2:Nf+1); yf(2:Nf+1) = yf(2:Nf+1) + dt*vf(2:Nf+1);  % move front
 xf(1) = xf(Nf+1); yf(1) = yf(Nf+1); xf(Nf+2) = xf(2); yf(Nf+2) = yf(2); 
 	
 d(2:nx+1,2:ny+1)=4*Lx;                                        % placeholder distance to front
 for l=2:Nf+1
   nfx = -(yf(l+1)-yf(l)); nfy =  (xf(l+1)-xf(l));             % outer normal vector
   ds  = sqrt(nfx*nfx+nfy*nfy); nfx = nfx/ds; nfy = nfy/ds;    % unit outer normal
   xfront = 0.5*(xf(l)+xf(l+1)); yfront = 0.5*(yf(l)+yf(l+1)); % find midpoint in the segment
   ip=floor((xfront+0.5*dx)/dx)+1; jp=floor((yfront+0.5*dy)/dy)+1;

   d1=sqrt(((xfront-  x(ip)))^2+((yfront-  y(jp)))^2);
   d2=sqrt(((xfront-x(ip+1)))^2+((yfront-  y(jp)))^2);
   d3=sqrt(((xfront-x(ip+1)))^2+((yfront-y(jp+1)))^2);
   d4=sqrt(((xfront-  x(ip)))^2+((yfront-y(jp+1)))^2);

   if d1<d(ip,jp)
    d(ip,jp) = d1; dn1 = (x(ip) - xfront)*nfx/dx + (y(jp) - yfront)*nfy/dy;
    chi(ip,jp) = 0.5*(1 + sign(dn1)); 
    if abs(dn1) < 0.5; chi(ip,jp) = 0.5 + dn1; end;
   end
   if d2<d(ip+1,jp)
    d(ip+1,jp) = d2; dn2 = (x(ip+1) - xfront)*nfx/dx + (y(jp) - yfront)*nfy/dy;
    chi(ip+1,jp) = 0.5*(1 + sign(dn2));
    if abs(dn2) < 0.5; chi(ip+1,jp) = 0.5 + dn2; end;
   end
   if d3<d(ip+1,jp+1)
    d(ip+1,jp+1) = d3; dn3 = (x(ip+1)-xfront)*nfx/dx+(y(jp+1)-yfront)*nfy/dy;
    chi(ip+1,jp+1)=0.5*(1 + sign(dn3)); 
    if abs(dn3) < 0.5; chi(ip+1,jp+1) = 0.5 + dn3; end;
   end
   if d4<d(ip,jp+1)
    d(ip,jp+1) = d4; dn4 = (x(ip) - xfront)*nfx/dx + (y(jp+1) - yfront)*nfy/dy;
    chi(ip,jp+1) = 0.5*(1 + sign(dn4)); 
    if abs(dn4) < 0.5; chi(ip,jp+1) = 0.5 + dn4; end;
   end
  end
              
  ro = r;
  r  = rho1*chi + rho2*(1-chi);  % obtain density from charact func
  mu = mu1*chi + mu2*(1-chi);
  cond = k1*chi + k2*(1-chi);
  rcp = rcp1*chi + rcp2*(1-chi);

  for l=1:Nf+1
    ds=sqrt((xf(l+1)-xf(l))^2+(yf(l+1)-yf(l))^2);
    tx(l)=(xf(l+1)-xf(l))/ds; ty(l)=(yf(l+1)-yf(l))/ds; % unit tangent vectors
  end
  tx(Nf+2)=tx(2); ty(Nf+2)=ty(2);

  fgx=zeros(nx+2,ny+2); fgy=zeros(nx+2,ny+2);
  for l=2:Nf+1                                               % distribute to the fixed grid
    fglx = gamma*(tx(l)-tx(l-1)); fgly = gamma*(ty(l)-ty(l-1));
    ip   = floor(xf(l)/dx)+1;     jp   = floor((yf(l)+0.5*dy)/dy)+1;
    ax   = xf(l)/dx-ip+1;         ay   = (yf(l)+0.5*dy)/dy-jp+1;
    fgx(ip,jp)     = fgx(ip,jp)     + (1.0-ax)*(1.0-ay)*fglx/dx/dy;
    fgx(ip+1,jp)   = fgx(ip+1,jp)   + ax*(1.0-ay)*fglx/dx/dy;
    fgx(ip,jp+1)   = fgx(ip,jp+1)   + (1.0-ax)*ay*fglx/dx/dy;
    fgx(ip+1,jp+1) = fgx(ip+1,jp+1) + ax*ay*fglx/dx/dy;

    ip = floor((xf(l)+0.5*dx)/dx)+1; jp = floor(yf(l)/dy)+1;
    ax = (xf(l)+0.5*dx)/dx-ip+1;     ay = yf(l)/dy-jp+1;	  
    fgy(ip,jp)     = fgy(ip,jp)     + (1.0-ax)*(1.0-ay)*fgly/dx/dy;
    fgy(ip+1,jp)   = fgy(ip+1,jp)   + ax*(1.0-ay)*fgly/dx/dy;
    fgy(ip,jp+1)   = fgy(ip,jp+1)   + (1.0-ax)*ay*fgly/dx/dy;
    fgy(ip+1,jp+1) = fgy(ip+1,jp+1) + ax*ay*fgly/dx/dy;  
  end

  fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
  fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  
  
  %u(1:nx+1,1) = u(1:nx+1,2); u(1:nx+1,ny+2) = u(1:nx+1,ny+1); % tangential vel BC
  %v(1,1:ny+1) =v(2,1:ny+1); v(nx+2,1:ny+1) = v(nx+1,1:ny+1); % tangential vel BC
  %v(1:nx+2,ny+1)= v(1:nx+2,ny); v(1:nx+2,1)= v(1:nx+2,2);
  %u(1,1:ny+2)= u(2,1:ny+2); u(nx+1,1:ny+2)= u(nx,1:ny+2); 
  u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
  v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC

  for i=2:nx; for j=2:ny+1     % temporary u-velocity (boundary values are not touched)
    ut(i,j) = (2.0/(r(i+1,j)+r(i,j)))*(0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+ dt* (...
            - (0.25/dx)*(ro(i+1,j)*(u(i+1,j)+u(i,j))^2-ro(i,j)*(u(i,j)+u(i-1,j))^2)...
            - (0.0625/dy)*( (ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1))*(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) ...
            - (ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))*(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))...
            + (2*mu(i+1,j)*(u(i+1,j)-u(i,j))/dx - 2*mu(i,j)*(u(i,j)-u(i-1,j))/dx)/dx ... 
            + (0.25*(mu(i,j)+mu(i+1,j)+mu(i,j+1)+mu(i+1,j+1))*((u(i,j+1)-u(i,j))/dy+ (v(i+1,j)-v(i,j))/dx) ...
            -0.25*(mu(i,j)+mu(i+1,j)+mu(i,j-1)+mu(i+1,j-1))*((u(i,j)-u(i,j-1))/dy+ (v(i+1,j-1)-v(i,j-1))/dx))/dy...
            + 0.5*(ro(i+1,j)+ro(i,j))*gx + fgx(i,j) ) );
  end; end

  for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
    vt(i,j) = (2.0/(r(i,j+1)+r(i,j)))*(0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+ dt* (...     
            - (0.0625/dx)*( (ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) ...
            - (ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))*(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
            - (0.25/dy)*(ro(i,j+1)*(v(i,j+1)+v(i,j))^2-ro(i,j)*(v(i,j)+v(i,j-1))^2 )...
             + (2*mu(i,j+1)*(v(i,j+1)-v(i,j))/dy - 2*mu(i,j)*(v(i,j)-v(i,j-1))/dy)/dy ... 
            + (0.25*(mu(i,j)+mu(i+1,j)+mu(i,j+1)+mu(i+1,j+1))*((u(i,j+1)-u(i,j))/dy+ (v(i+1,j)-v(i,j))/dx) ...
            -0.25*(mu(i,j)+mu(i,j+1)+mu(i-1,j+1)+mu(i-1,j))*((u(i-1,j+1)-u(i-1,j))/dy+ (v(i,j)-v(i-1,j))/dx))/dy...
            + 0.5*(ro(i,j+1)+ro(i,j))*gy + fgy(i,j) ) );    
  end; end     

  for i = 2:nx+1; for j = 2:ny+1
    tmp1(i,j) = (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy - 1/hlg*(1/rho2-1/rho1)*qdot(i,j));
    tmp2(i,j) =1/( (1/dx)*(1/(dx*(r(i+1,j)+r(i,j)))+ 1/(dx*(r(i-1,j)+r(i,j))) )+ ...
                   (1/dy)*(1/(dy*(r(i,j+1)+r(i,j)))+ 1/(dy*(r(i,j-1)+r(i,j))) )   );
  end; end

  for it = 1:maxit	               % solve for pressure by SOR
    pold   = p;
    p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
    for i=2:nx+1; for j=2:ny+1
      p(i,j) = (1.0-omg)*p(i,j) + omg*tmp2(i,j)*(        ...
      (1/dx)*( p(i+1,j)/(dx*(r(i+1,j)+r(i,j)))+ p(i-1,j)/(dx*(r(i-1,j)+r(i,j))) )+    ...
      (1/dy)*( p(i,j+1)/(dy*(r(i,j+1)+r(i,j)))+ p(i,j-1)/(dy*(r(i,j-1)+r(i,j))) ) - tmp1(i,j));
    end; end
    if max(max(abs(pold-p))) < maxError; break;  end
  end
                                      
  for i=2:nx; for j=2:ny+1   % correct the u-velocity 
    u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
  end; end
      
  for i=2:nx+1; for j=2:ny   % correct the v-velocity 
    v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
  end; end
  
  c(1,:)=c(2,:);c(nx+2,:)=c(nx+1,:);c(:,1)=c(:,2);c(:,ny+2)=c(:,ny+1);
  for i=2:nx+1; for j=2:ny+1
          c(i,j) = c(i,j)+dt*(-(u(i,j)*(c(i,j)+c(i+1,j))/2.0-u(i-1,j)*(c(i,j)+c(i-1,j))/2.0)/dx - ...
              (v(i,j)*(c(i,j)+c(i,j+1))/2.0-v(i,j-1)*(c(i,j)+c(i,j-1))/2.0)/dy + ...
              0.5/rcp(i,j)*( ( (cond(i+1,j)+cond(i,j)) * (c(i+1,j)-c(i,j))/dx - (cond(i-1,j)+cond(i,j)) * (c(i,j)-c(i-1,j))/dx ) /dx + ...
              ( (cond(i,j)+cond(i,j+1)) * (c(i,j+1)-c(i,j))/dy - (cond(i,j)+cond(i,j-1)) * (c(i,j)-c(i,j-1))/dy ) /dy ) ...
              - qdot(i,j)/rcp(i,j));
              
      end;end;
  xfold=xf; yfold=yf; j=1;   % add or/and delete points in the front
  for l=2:Nf+1
    ds=sqrt( ((xfold(l)-xf(j))/dx)^2 + ((yfold(l)-yf(j))/dy)^2);
    if (ds > 0.5)                                    % add point
      j=j+1; xf(j) = 0.5*(xfold(l)+xf(j-1)); yf(j) = 0.5*(yfold(l)+yf(j-1));
      j=j+1; xf(j) = xfold(l); yf(j) = yfold(l);
    elseif (ds < 0.25)                               % skip point
    else
      j = j+1; xf(j) = xfold(l); yf(j) = yfold(l);   % copy point    
    end    
  end
  Nf = j-1; xf(1) = xf(Nf+1); yf(1) = yf(Nf+1); xf(Nf+2) = xf(2); yf(Nf+2) = yf(2);

  
  time = time+dt
  
  if (mod(is,plot_freq)==0) | (is==1);                         % plot solution
  uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
  vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
  plotc = c;
  %frand=find(plotc(:)>Tw);
  %frand2=find(plotc(:)<Tin);
  %plotc(frand)=Tw;
  %plotc(frand2)=Tin;
%   fighandle=figure(1); 
%   contourf(x,y,plotc','linestyle','none'); colorbar; colormap(parula(500)); ax = gca; %set(ax,'CLim',[Tsat Tw]); 
%   hold on;
%    
%   h1=quiver(xh(1:2:end,1:2:end),yh(1:2:end,1:2:end),uu(1:2:end,1:2:end)',vv(1:2:end,1:2:end)');axis equal; 
%   axis([0 Lx 0 Ly]);hold on; 
%   set(h1,'AutoScale','on', 'AutoScaleFactor', 2);
%   hold on; plot(xf(1:Nf),yf(1:Nf),'k','linewidth',3); drawnow; hold off
  %filename = ['C:\Users\ganeshm\OneDrive - purdue.edu\lenovo_backup\ME697\figures_hw3\' sprintf('P1_%d',is) '.png'];
  %uncomment to save
  %saveas(fighandle,filename);
  %filename = ['C:\Users\ganeshm\OneDrive - purdue.edu\lenovo_backup\ME697\figures_hw3\' sprintf('P1_%d',is) '.eps'];
  %saveas(fighandle,filename);
  end

end   
    end
end

%Plotting for Tw variation
timed=linspace(0,nstep*dt,nstep);
figure;plot(timed,squeeze(d_n(:,1:length(Twvec),1)),'linewidth',3)
legend('(T_\infty-T_{sat})/T_{sat}=9','(T_\infty-T_{sat})/T_{sat}=4','(T_\infty-T_{sat})/T_{sat}=1','location','northwest')
xlabel('Time')
ylabel('Diameter')

%uncomment for Plotting for kl variation
% figure;plot(timed,squeeze(d_n(:,1,1:length(klvec))),'linewidth',3)
% legend('k_l/k_g=1','k_l/k_g=2.5','k_l/k_g=5.0','location','northwest')
% xlabel('Time')
% ylabel('Diameter')

%Plotting mu contour
fighandle=figure(1); 
contourf(x,y,mu','linestyle','none'); colorbar; colormap(jet(500)); ax = gca; %set(ax,'CLim',[Tsat Tw]); 
hold on;
%    
h1=quiver(xh(1:2:end,1:2:end),yh(1:2:end,1:2:end),uu(1:2:end,1:2:end)',vv(1:2:end,1:2:end)','color','#2F4F4F');axis equal; 
axis([0 Lx 0 Ly]);hold on; 
set(h1,'AutoScale','on', 'AutoScaleFactor', 2);
hold on; plot(xf(1:Nf),yf(1:Nf),'k','linewidth',3); drawnow; hold off
%filename = ['C:\Users\ganeshm\OneDrive - purdue.edu\lenovo_backup\ME697\project\variable_q\' sprintf('mu_var.png')];
%print(filename,'-dpng','-r400')
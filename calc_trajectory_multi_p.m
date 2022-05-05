function [floes,allF]=calc_trajectory_multi_p(dt,ocean,winds,floes, floeShapes, activeFloes)

ext_force=[0 0];
ext_torque=0;

Xo=ocean.Xo;
Yo=ocean.Yo;

Xw = winds.Xw*1e3;
Yw = winds.Yw*1e3;

XoMax = max(max(Xo));
XoMin = min(min(Xo));
YoMax = max(max(Yo));
YoMin = min(min(Yo));

XwMax = max(max(Xw));
XwMin = min(min(Xw));
YwMax = max(max(Yw));
YwMin = min(min(Yw));

Uocn=ocean.Uocn;
Vocn=ocean.Vocn;

Uw = winds.Uw;
Vw = winds.Vw;
%[Xocn, Yocn]=meshgrid(Xo,Yo);
dXo=Xo(2)-Xo(1);

nFloes = length(floes);

for iFloe = 1:nFloes
    if ~activeFloes(iFloe)
        continue
    end
    
    floeXi = floes(iFloe).Xi;
    floeYi = floes(iFloe).Yi;
    floeUi = floes(iFloe).Ui;
    floeVi = floes(iFloe).Vi;
    
    Xi=floeXi;
    Yi=floeYi;
     


    % Spin up of an ice floe by an ocean eddy;
    % ice floe is a disk with a radius R_floe;
    % ocean eddy has a uniform vorticity and is
    % centered at the center of the ice floe (hence no drift for now).

    % ice-water drag is parameterized as tau=rho0*Cd*|Uocean-Uice|*(Uocean-Uice)
    % no turning angle is used for now

    % ice-ocean parameters

    rho_ice=920; % kg/m3
    rho0=1027;   % ocean density
    Cd=5.5e-3;
    
    % ice-water drag coefficient
    rho_air=1.2;
    Cd_atm=1.6e-3;

    fc=1.4e-4; %coriolis parameter

    floe_area=floes(iFloe).area;
    floe_mass=floes(iFloe).mass; % total mass
    floe_thickness=floe_mass/floes(iFloe).area/rho_ice;
    floe_inertia_moment=floes(iFloe).inertia_moment; % moment of inertia
    R_floe=floes(iFloe).rmax;

    A_alpha = squeeze(floeShapes{iFloe}(mod(fix(-floes(iFloe).alpha_i/pi*180), 360)+1, :, :));
    floe_mask=(A_alpha==1);

    Xg=floes(iFloe).X+floeXi; % grid centered around the ice floe
    Yg=floes(iFloe).Y+floeYi;

    [theta,rho] = cart2pol(Xg-Xi,Yg-Yi);

    Uice = floeUi;
    Vice = floeVi;

    Uocn_interp=interp2(Xo,Yo, Uocn, mod(Xg-XoMin, XoMax-XoMin)+XoMin, mod(Yg-YoMin, YoMax-YoMin)+YoMin, 'linear', 0);
    Vocn_interp=interp2(Xo,Yo, Vocn, mod(Xg-XoMin, XoMax-XoMin)+XoMin, mod(Yg-YoMin, YoMax-YoMin)+YoMin, 'linear', 0);

    Fx_pressureGrad=-rho_ice*floe_thickness*fc*Vocn_interp;
    Fy_pressureGrad=+rho_ice*floe_thickness*fc*Uocn_interp;

    Uw_interp=interp2(Xw,Yw, Uw, mod(Xi-XoMin, XoMax-XoMin)+XoMin, mod(Yi-YoMin, YoMax-YoMin)+YoMin, 'linear', 0);
    Vw_interp=interp2(Xw,Yw, Vw, mod(Xi-XoMin, XoMax-XoMin)+XoMin, mod(Yi-YoMin, YoMax-YoMin)+YoMin, 'linear', 0);
    
    duw=Uw_interp-Uice; 
    dvw=Vw_interp-Vice;
    
    turn_angle_atm=0*pi/180;

    Fx_atm=rho_air*Cd_atm*sqrt(duw.^2+dvw.^2).*( cos(turn_angle_atm)*duw+sin(turn_angle_atm)*dvw);
    Fy_atm=rho_air*Cd_atm*sqrt(duw.^2+dvw.^2).*(-sin(turn_angle_atm)*duw+cos(turn_angle_atm)*dvw);

    du=Uocn_interp-Uice; 
    dv=Vocn_interp-Vice;

    turn_angle_ocn=20*pi/180;

    tau_ocnX=rho0*Cd*sqrt(du.^2+dv.^2).*( cos(turn_angle_ocn)*du+sin(turn_angle_ocn)*dv);
    tau_ocnY=rho0*Cd*sqrt(du.^2+dv.^2).*(-sin(turn_angle_ocn)*du+cos(turn_angle_ocn)*dv);

    Fx=tau_ocnX+Fx_atm+Fx_pressureGrad; % add up forces
    Fy=tau_ocnY+Fy_atm+Fy_pressureGrad;

    % updating the ice floe vorticity with averaged torques over the ice floe area
    torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque

    %adding coriolis force; it has no torque.
    Fx=Fx+rho_ice*floe_thickness*fc*floeVi;
    Fy=Fy-rho_ice*floe_thickness*fc*floeUi;

    % updating the ice floe coordinates with velocities
    floeXi=floeXi+dt*floeUi;
    floes(iFloe).dXi_p=floeUi;

    floeYi=floeYi+dt*floeVi;
    floes(iFloe).dYi_p=floeVi;
    
    floes(iFloe).alpha_i=floes(iFloe).alpha_i+dt*floes(iFloe).ksi_ice;
    floes(iFloe).dalpha_i_p=floes(iFloe).ksi_ice;
    

    dUi_dt=(mean(Fx(floe_mask))*floe_area+ext_force(1))/floe_mass;
    floeUi=floeUi+dt*dUi_dt;  floes(iFloe).dUi_p=dUi_dt;

    dVi_dt=(mean(Fy(floe_mask))*floe_area+ext_force(2))/floe_mass;
    floeVi=floeVi+dt*dVi_dt;  floes(iFloe).dVi_p=dVi_dt;

    dksi_ice_dt=(mean(torque(floe_mask))*floe_area+ext_torque)/floe_inertia_moment;
    floes(iFloe).ksi_ice=floes(iFloe).ksi_ice+dt*dksi_ice_dt;
    floes(iFloe).dksi_ice_p=dksi_ice_dt;

    floes(iFloe).Xi = floeXi;
    floes(iFloe).Yi = floeYi;
    floes(iFloe).Ui = floeUi;
    floes(iFloe).Vi = floeVi;

end

end



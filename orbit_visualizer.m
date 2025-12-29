function orbit_visualizer()
% ORBIT_VISUALIZER Interactive 3D Elliptical Orbit Visualizer (Earth centered, ECI frame)
%
% Drag sliders to change orbital elements in real time:
% a (km), e, i (deg), RAAN Ω (deg), ω argument of periapsis (deg), ν true anomaly (deg)
% Shows full orbit path, current spacecraft position, and Earth as a sphere
%
% Run: orbit_visualizer
%
% Author: Grace Fermelia

%% Constants
R_EARTH = 6378.137;      % km
NUM_POINTS = 800;

%% Default orbital elements
a0    = max(12000, R_EARTH + 150);   % semi major axis
e0    = max(0.0, min(0.95, 0.3));    % eccentricity clamped
i0    = deg2rad(30.0);               % inclination
raan0 = deg2rad(40.0);               % RAAN
argp0 = deg2rad(60.0);               % argument of periapsis
nu0   = deg2rad(0.0);                % true anomaly

%% Figure and axes
fig = figure('Name','Interactive 3D Orbit Visualizer','NumberTitle','off',...
    'Position',[100 100 1200 800]);

ax = axes('Position',[0.25 0.18 0.70 0.72]);
hold(ax,'on');
grid(ax,'on');
axis(ax,'equal');
view(ax,3);

%% Initial orbit
nus_full = linspace(0,2*pi,NUM_POINTS);
[x,y,z] = orbit_xyz(a0,e0,i0,raan0,argp0,nus_full);

make_earth(ax);
make_eci_axes(ax);

orbit_line = plot3(ax,x,y,z,'b','LineWidth',2);

[x0,y0,z0] = orbit_xyz(a0,e0,i0,raan0,argp0,nu0);
sc_marker = scatter3(ax,x0,y0,z0,60,'filled','MarkerFaceColor','red');

[node_start,node_end,apses_start,apses_end] = make_orbital_lines(a0,e0,i0,raan0,argp0);

line_of_nodes = plot3(ax,[node_start(1) node_end(1)],...
    [node_start(2) node_end(2)],...
    [node_start(3) node_end(3)],'--','LineWidth',2,'Color',[1 0.5 0]);

line_of_apses = plot3(ax,[apses_start(1) apses_end(1)],...
    [apses_start(2) apses_end(2)],...
    [apses_start(3) apses_end(3)],'--','LineWidth',2,'Color',[0.5 0 0.5]);

radius_vector = quiver3(ax,0,0,0,x0,y0,z0,0,...
    'Color','cyan','LineWidth',3,'MaxHeadSize',0.05);

xlabel(ax,'X (km)');
ylabel(ax,'Y (km)');
zlabel(ax,'Z (km)');
title(ax,format_title(a0,e0,i0,raan0,argp0),'Interpreter','none');

set_equal_3d(ax,[x;y;z]);

%% Sliders
s_a = create_vertical_slider([0.05 0.70 0.02 0.20],'a (km)',R_EARTH*1.02,60000,a0);
s_e = create_vertical_slider([0.05 0.45 0.02 0.20],'e',0,0.95,e0);
s_i = create_vertical_slider([0.05 0.20 0.02 0.20],'i (deg)',0,180,rad2deg(i0));

s_raan = create_horizontal_slider([0.30 0.05 0.25 0.03],'Ω (deg)',0,360,rad2deg(raan0));
s_argp = create_horizontal_slider([0.60 0.05 0.25 0.03],'ω (deg)',0,360,rad2deg(argp0));
s_nu   = create_horizontal_slider([0.30 0.11 0.55 0.03],'ν (deg)',0,360,rad2deg(nu0));

uicontrol('Style','pushbutton','String','Reset',...
    'Units','normalized','Position',[0.05 0.05 0.10 0.05],...
    'Callback',@reset_callback);

%% Store handles
handles = struct();
handles.ax = ax;
handles.orbit_line = orbit_line;
handles.sc_marker = sc_marker;
handles.line_of_nodes = line_of_nodes;
handles.line_of_apses = line_of_apses;
handles.radius_vector = radius_vector;
handles.sliders = struct('a',s_a,'e',s_e,'i',s_i,'raan',s_raan,'argp',s_argp,'nu',s_nu);
handles.nus_full = nus_full;
set(fig,'UserData',handles);

set([s_a s_e s_i s_raan s_argp s_nu],'Callback',@(s,e)update_orbit(fig));

end

%% Callback functions
function update_orbit(fig)
handles = get(fig,'UserData');

a    = get(handles.sliders.a,'Value');
e    = get(handles.sliders.e,'Value');
i    = deg2rad(get(handles.sliders.i,'Value'));
raan = deg2rad(get(handles.sliders.raan,'Value'));
argp = deg2rad(get(handles.sliders.argp,'Value'));
nu   = deg2rad(get(handles.sliders.nu,'Value'));

[x,y,z] = orbit_xyz(a,e,i,raan,argp,handles.nus_full);
set(handles.orbit_line,'XData',x,'YData',y,'ZData',z);

[xm,ym,zm] = orbit_xyz(a,e,i,raan,argp,nu);
set(handles.sc_marker,'XData',xm,'YData',ym,'ZData',zm);

delete(handles.radius_vector);
handles.radius_vector = quiver3(handles.ax,0,0,0,xm,ym,zm,0,...
    'Color','cyan','LineWidth',3,'MaxHeadSize',0.05);

[node_start,node_end,apses_start,apses_end] = make_orbital_lines(a,e,i,raan,argp);

set(handles.line_of_nodes,'XData',[node_start(1) node_end(1)],...
    'YData',[node_start(2) node_end(2)],...
    'ZData',[node_start(3) node_end(3)]);

set(handles.line_of_apses,'XData',[apses_start(1) apses_end(1)],...
    'YData',[apses_start(2) apses_end(2)],...
    'ZData',[apses_start(3) apses_end(3)]);

set_equal_3d(handles.ax,[x;y;z]);
title(handles.ax,format_title(a,e,i,raan,argp),'Interpreter','none');

set(fig,'UserData',handles);
end

function reset_callback(~,~)
fig = gcbf;
handles = get(fig,'UserData');

set(handles.sliders.a,'Value',12000);
set(handles.sliders.e,'Value',0.3);
set(handles.sliders.i,'Value',30);
set(handles.sliders.raan,'Value',40);
set(handles.sliders.argp,'Value',60);
set(handles.sliders.nu,'Value',0);

update_orbit(fig);
end

%% Math functions
function [x,y,z] = orbit_xyz(a,e,i,raan,argp,nus)
p = a*(1-e^2);
r = p./(1+e*cos(nus));

x_p = r.*cos(nus);
y_p = r.*sin(nus);
z_p = zeros(size(nus));

Q = pqw_to_eci_matrix(raan,i,argp);
r_eci = Q*[x_p;y_p;z_p];

x = r_eci(1,:);
y = r_eci(2,:);
z = r_eci(3,:);
end

function Q = pqw_to_eci_matrix(raan,inc,argp)
R3_raan = [cos(raan) -sin(raan) 0; sin(raan) cos(raan) 0; 0 0 1];
R1_inc  = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
R3_argp = [cos(argp) -sin(argp) 0; sin(argp) cos(argp) 0; 0 0 1];
Q = R3_raan*R1_inc*R3_argp;
end

function make_earth(ax)
R = 6378.137;
[xs,ys,zs] = sphere(32);
surf(ax,xs*R,ys*R,zs*R,'FaceAlpha',0.15,'EdgeColor','none','FaceColor','blue');
end

function make_eci_axes(ax)
R = 6378.137*2;
quiver3(ax,0,0,0,R,0,0,0,'r','LineWidth',2);
quiver3(ax,0,0,0,0,R,0,0,'g','LineWidth',2);
quiver3(ax,0,0,0,0,0,R,0,'b','LineWidth',2);
end

function [ns,ne,as,ae] = make_orbital_lines(a,e,i,raan,argp)
L = 2*a;
ns = -L*[cos(raan);sin(raan);0];
ne =  L*[cos(raan);sin(raan);0];
Q = pqw_to_eci_matrix(raan,i,argp);
d = Q*[1;0;0];
as = -L*d;
ae =  L*d;
end

function set_equal_3d(ax,xyz)
x=xyz(1,:); y=xyz(2,:); z=xyz(3,:);
r = max([range(x) range(y) range(z)])/2;
m = [mean(x) mean(y) mean(z)];
xlim(ax,m(1)+[-r r]); ylim(ax,m(2)+[-r r]); zlim(ax,m(3)+[-r r]);
end

function t = format_title(a,e,i,raan,argp)
t = sprintf('a=%.0f km, e=%.2f, i=%.1f°, Ω=%.1f°, ω=%.1f°',...
    a,e,rad2deg(i),rad2deg(raan),rad2deg(argp));
end

function s = create_vertical_slider(pos,label,minv,maxv,val)
s = uicontrol('Style','slider','Units','normalized','Position',pos,...
    'Min',minv,'Max',maxv,'Value',val);
uicontrol('Style','text','Units','normalized',...
    'Position',[pos(1)-0.01 pos(2)+pos(4)+0.01 0.05 0.03],...
    'String',label,'FontSize',8);
end

function s = create_horizontal_slider(pos,label,minv,maxv,val)
s = uicontrol('Style','slider','Units','normalized','Position',pos,...
    'Min',minv,'Max',maxv,'Value',val);
uicontrol('Style','text','Units','normalized',...
    'Position',[pos(1) pos(2)-0.02 pos(3) 0.02],...
    'String',label,'FontSize',8);
end

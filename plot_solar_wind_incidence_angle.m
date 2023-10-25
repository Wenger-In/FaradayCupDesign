clear; close all;
%%
file_dir = 'E:\Research\Data\PSP\';
file_name = 'SPP_(HCI)_(20250101T000000-20250831T000000)_copy.dat';
pos_file = [file_dir, file_name];
pos_data = importdata(pos_file);
pos_data = pos_data(500:850,:); % select time range
[num_epoch,~] = size(pos_data);

index = 1 : num_epoch;
year = pos_data(:,1);
month = pos_data(:,2);
day = pos_data(:,3);
hour = pos_data(:,4);
minute = pos_data(:,5);
second = pos_data(:,6);
microsecond = pos_data(:,7);
sc_pos_HCIx = pos_data(:,8); % [km]
sc_pos_HCIy = pos_data(:,9);
sc_pos_HCIz = pos_data(:,10);
sc_vel_HCIx = pos_data(:,11); % [km/s]
sc_vel_HCIy = pos_data(:,12);
sc_vel_HCIz = pos_data(:,13);
%%
[sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn] = calc_HCI2SCRTN(sc_vel_HCIx,sc_vel_HCIy,sc_vel_HCIz,sc_pos_HCIx,sc_pos_HCIy,sc_pos_HCIz);
sc_pos_mod = sqrt(sc_pos_HCIx.^2 + sc_pos_HCIy.^2 + sc_pos_HCIz.^2);
%%
Rs = 6.955e5; % [km]
V_sw = 100; % [km/s]
V_sw2sc_RTNr = V_sw - sc_vel_RTNr;
V_sw2sc_RTNt = -sc_vel_RTNt;
V_sw2sc_RTNn = -sc_vel_RTNn;
V_sw2sc_mod = sqrt(V_sw2sc_RTNr.^2 + V_sw2sc_RTNt.^2 + V_sw2sc_RTNn.^2);
%%
angle = zeros(num_epoch,1);
for i_epoch = 1 : num_epoch
    angle = acosd(V_sw2sc_RTNr./V_sw2sc_mod);
end
%%
figure();

subplot(4,1,1)
% plot(sc_pos_mod/Rs,V_sw2sc_RTNr,'LineWidth',2); grid on
patch(sc_pos_mod/Rs,V_sw2sc_RTNr,index,index,'LineWidth',4,'Edgecolor','flat','Facecolor','none'); grid on
c = colorbar; colormap jet
c.Ticks = [125,175,225];
c.TickLabels = {'06/06','06/18','06/30'};
% xlabel('r [Rs]');
ylabel('V_{sw2sc,R} [km]');
set(gca,'LineWidth',2,'FontSize',15,'XTickLabel','','Clim',[125,225]);

subplot(4,1,2)
% plot(sc_pos_mod/Rs,V_sw2sc_RTNt,'LineWidth',2); grid on
patch(sc_pos_mod/Rs,V_sw2sc_RTNt,index,index,'LineWidth',4,'Edgecolor','flat','Facecolor','none'); grid on
c = colorbar; colormap jet
c.Ticks = [125,175,225];
c.TickLabels = {'06/06','06/18','06/30'};
% xlabel('r [Rs]');
ylabel('V_{sw2sc,T} [km]');
set(gca,'LineWidth',2,'FontSize',15,'XTickLabel','','Clim',[125,225]);

subplot(4,1,3)
% plot(sc_pos_mod/Rs,V_sw2sc_RTNn,'LineWidth',2); grid on
patch(sc_pos_mod/Rs,V_sw2sc_RTNn,index,index,'LineWidth',4,'Edgecolor','flat','Facecolor','none'); grid on
c = colorbar; colormap jet
c.Ticks = [125,175,225];
c.TickLabels = {'06/06','06/18','06/30'};
% xlabel('r [Rs]');
ylabel('V_{sw2sc,N} [km]');
set(gca,'LineWidth',2,'FontSize',15,'XTickLabel','','Clim',[125,225]);

subplot(4,1,4)
% plot(sc_pos_mod/Rs,angle,'LineWidth',2); grid on
patch(sc_pos_mod/Rs,angle,index,index,'LineWidth',4,'Edgecolor','flat','Facecolor','none'); grid on
c = colorbar; colormap jet
c.Ticks = [125,175,225];
c.TickLabels = {'06/06','06/18','06/30'};
xlabel('r [Rs]');
ylabel('\theta_{i} [deg.]');
set(gca,'LineWidth',2,'FontSize',15,'Clim',[125,225]);

sgtitle(['20250506-20250801 V\_sw=',num2str(V_sw),'km/s'],'FontSize',20);
%%

%% functions
function [x_RTN,y_RTN,z_RTN] = calc_HCI2SCRTN(x_HCI,y_HCI,z_HCI,SC_HCIx,SC_HCIy,SC_HCIz)
% Change the coordiantes xyz in HCI frame to xyz in spacecraft RTN frame
%   input: x_HCI, y_HCI, z_HCI, the velocity in HCI frame (km/s)
%          SC_HCIx, SC_HCIy, SC_HCIz, the sapcecraft position in HCI frame (km)
%   output: x_RTN, y_RTN, z_RTN,the velocity in SC RTN frame (km/s)
%   This function does not consider the move of the origin (using for velocity conversion)
    num = length(x_HCI);
    xyz_RTN = zeros(num,3);
    for i = 1:num
        Q = zeros(3,3);
        x1 = [1 0 0];
        y1 = [0 1 0];
        z1 = [0 0 1];
        x2 = [SC_HCIx(i),SC_HCIy(i),SC_HCIz(i)];
        if norm(x2)~= 0
            x2 = x2/norm(x2);
        end
        y2 = cross(z1,x2);
        if norm(y2)~= 0
            y2 = y2/norm(y2);
        end
        z2 = cross(x2,y2);
        if norm(z2)~= 0
            z2 = z2/norm(z2);
        end
        Q(1,1) = dot(x2,x1); Q(1,2) = dot(x2,y1); Q(1,3) = dot(x2,z1);
        Q(2,1) = dot(y2,x1); Q(2,2) = dot(y2,y1); Q(2,3) = dot(y2,z1);
        Q(3,1) = dot(z2,x1); Q(3,2) = dot(z2,y1); Q(3,3) = dot(z2,z1);
        xyz_RTN(i,:) = Q*[x_HCI(i);y_HCI(i);z_HCI(i)];
    end
    x_RTN = xyz_RTN(:,1); y_RTN = xyz_RTN(:,2); z_RTN = xyz_RTN(:,3);
end
clear; close all;
data_store = 'E:\Research\Data\PSP\Encounter 7\';
DataDir = dir([data_store 'psp_swp_spc_*.cdf']);
%%
Rs = 6.955e8;
sc_pos_mod_tot = [];
sc_vel_RTNr_tot = [];
sc_vel_RTNt_tot = [];
sc_vel_RTNn_tot = [];
for i = 1 : length(DataDir)
    spc_file = [data_store,DataDir(i).name];
    spc_info = spdfcdfinfo(spc_file);
    %%
    epoch = spdfcdfread(spc_file,'Variables','Epoch');
    general_flag = spdfcdfread(spc_file,'Variables','general_flag');
    sc_pos_HCI = spdfcdfread(spc_file,'Variables','sc_pos_HCI');
    sc_vel_HCI = spdfcdfread(spc_file,'Variables','sc_vel_HCI');
    np_moment = spdfcdfread(spc_file,'Variables','np_moment');
    %%
    sc_pos_HCIx = sc_pos_HCI(:,1); sc_pos_HCIy = sc_pos_HCI(:,2); sc_pos_HCIz = sc_pos_HCI(:,3);
    sc_vel_HCIx = sc_vel_HCI(:,1); sc_vel_HCIy = sc_vel_HCI(:,2); sc_vel_HCIz = sc_vel_HCI(:,3);
    [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn] = calc_HCI2SCRTN(sc_vel_HCIx,sc_vel_HCIy,sc_vel_HCIz,sc_pos_HCIx,sc_pos_HCIy,sc_pos_HCIz);
    sc_pos_mod = sqrt(sc_pos_HCIx.^2 + sc_pos_HCIy.^2 + sc_pos_HCIz.^2);
%     sc_vel_mod = sqrt(sc_vel_HCIx.^2 + sc_vel_HCIy.^2 + sc_vel_HCIz.^2);
%     angle = zeros(length(epoch),1);
%     for j = 1 : length(epoch)
%         angle(j) = acosd(dot(sc_pos_HCI(j,:),sc_vel_HCI(j,:))/sc_pos_mod(j)/sc_vel_mod(j));
%     end
    %%
    sc_pos_mod_tot = [sc_pos_mod_tot; sc_pos_mod];
    sc_vel_RTNr_tot = [sc_vel_RTNr_tot; sc_vel_RTNr]; 
    sc_vel_RTNt_tot = [sc_vel_RTNt_tot; sc_vel_RTNt]; 
    sc_vel_RTNn_tot = [sc_vel_RTNn_tot; sc_vel_RTNn]; 
end
%%
V_sw = 200;
V_sw2sc_RTNr_tot = V_sw - sc_vel_RTNr_tot;
V_sw2sc_RTNt_tot = -sc_vel_RTNt_tot;
V_sw2sc_RTNn_tot = -sc_vel_RTNn_tot;
V_sw2sc_mod_tot = sqrt(V_sw2sc_RTNr_tot.^2 + V_sw2sc_RTNt_tot.^2 + V_sw2sc_RTNn_tot.^2);
%%
figure();
plot(sc_pos_mod_tot/Rs, V_sw2sc_mod_tot);
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
function [data_rotated] = scan2phys(data,ps,pm)
%% load an arbitrary set of points
normalmeas = cross(pm(2,:)-pm(1,:), pm(3,:)-pm(1,:)); % normal to the plane formed by the measured landmark coordinates
normalmeas = normalmeas/norm(normalmeas);
normalscan = cross(ps(2,:)-ps(1,:), ps(3,:)-ps(1,:)); % normal to the plane formed by the measured landmark coordinates
normalscan = normalscan/norm(normalscan);
vm = pm(2,:)-pm(1,:); % need to be orthogonal to normalmes
vm = vm/norm(vm);
vs = ps(2,:)-ps(1,:); % need to be orthogonal to normalscan
vs = vs/norm(vs);
%% define SCANNED coordinate system
scan_sys.nhat = normalscan; % \hat{n} = \hat{x}n_x+\hat{y}n_y+\hat{z}n_z
scan_sys.vhat = vs; % define scanned coordinate system orientation
scan_sys.p0 = ps(1,:); % p_0 = \hat{x}x_0+\hat{y}y_0+\hat{z}z_0 is the origin of the local coordinate system
[scan_sys.T,scan_sys.Ttilde] = PlaneCoordinatesRotationMatrices2(scan_sys.nhat,scan_sys.vhat); 
%% define MEASURED coordinate system
meas_sys.nhat = normalmeas; % \hat{n} = \hat{x}n_x+\hat{y}n_y+\hat{z}n_z
meas_sys.vhat = vm; % define measured coordinate system orientation
meas_sys.p0 = pm(1,:); % p_0 = \hat{x}x_0+\hat{y}y_0+\hat{z}z_0 is the origin of the local coordinate system
[meas_sys.T,meas_sys.Ttilde] = PlaneCoordinatesRotationMatrices2(meas_sys.nhat,meas_sys.vhat); 
%% transform SCANNED points to MEASURED coordinate system
fun_scan2meas = @(p) ((p-scan_sys.p0)*scan_sys.T)*meas_sys.Ttilde+meas_sys.p0;
data_rotated = fun_scan2meas(data);

end
% %% points in global coordinate system
% fn=1
% set(figure(fn),'color','white'); fn=fn+1; hold on;
% grid on; grid minor;
% view(3);
% daspect([1,1,1]);
% [nhat_plt,vhat_plt,pts_plt] = deal(cell(1,length(plane)));
% nhat_caption = [];
% vhat_caption = [];
% pts_caption = [];
% for i = 1:length(plane)
%     nhat_plt{i} = quiver3(plane{i}.p0(1),plane{i}.p0(2),plane{i}.p0(3),plane{i}.nhat(1),plane{i}.nhat(2),plane{i}.nhat(3),'-','linewidth',2,'MaxHeadSize',5);
%     vhat_plt{i} = quiver3(plane{i}.p0(1),plane{i}.p0(2),plane{i}.p0(3),plane{i}.vhat(1),plane{i}.vhat(2),plane{i}.vhat(3),'-','linewidth',2,'MaxHeadSize',5);
%     pts_plt{i} = plot3(points_global{i}(:,1),points_global{i}(:,2),points_global{i}(:,3),'x');
%     nhat_caption = [nhat_caption,{['Normal Plane ',num2str(i)]}];
%     vhat_caption = [vhat_caption,{['Orientation ',num2str(i)]}];
%     pts_caption = [pts_caption,{['Points Plane ',num2str(i)]}];
% end
% title('Set of points in global coordinates','interpreter','latex');
% xlabel('$x$ [m]','interpreter','latex');
% ylabel('$y$ [m]','interpreter','latex');
% zlabel('$z$ [m]','interpreter','latex');
% legend([[nhat_plt{:}],[vhat_plt{:}],[pts_plt{:}]],[nhat_caption,vhat_caption,pts_caption],'interpreter','latex','location','best');
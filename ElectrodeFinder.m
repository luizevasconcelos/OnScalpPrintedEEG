%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Authors: Sandhya Tiku, Luize Vasconcelos                                   %
% Nanshu Lu Research Group, UT austin                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [electrodepaths]=ElectrodeFinder(alltraces,nelectdisplay,diameter)
nelect=nelectdisplay;
% create flat path
nlines=6; %lines across diameter = nlines+1
R= diameter/2; %mm electrode radius
[xyzflat] = createflatpath(nlines,R);
electrodepaths=cell(nelect,1);

for j=1:nelect
    electrode(j,1:3)=alltraces{j}(1,1:3);
    normal(j,1:3)=alltraces{j}(1,4:6);
    [electrodepaths{j}] = moveElectrode(xyzflat,normal(j,:),electrode(j,:));
end

%electrodepaths cell contains each electrode path [X Y Z Nx Ny Nz]

    function path = createflatpath(nlines,R)
     
        %%option one
        resol=2*R/(nlines);
        x0=0;
        y0=-R;
        inc=(0:1:nlines)';
        y=inc*resol-R;
        y=[y;y];
        y=sort(y);
        for i=1:length(y)
            if mod(i,2)~=1 %odd
                x(i,1)=-sqrt(R.^2-y(i).^2);
            else
                x(i,1)=sqrt(R.^2-y(i).^2);
            end
        end
        xc=x(2:end-1);        
        yc=y(2:end-1);
        xc=[xc(3);xc;xc(end-2)];%close circle
        yc=[yc(2);yc;yc(end-1)];%close circle
    
         path= [xc yc zeros(length(yc),1)];

    end

    function points_global = moveElectrode(data,normal,electrode)
        
        %% load points defined in a local coordinate system
        points_local = data; % points in local coordinates
        %% define plane using normal and local origin
        nhat = normal; % \hat{n} = \hat{x}n_x+\hat{y}n_y+\hat{z}n_z
        p0 = electrode; % p_0 = \hat{x}x_0+\hat{y}y_0+\hat{z}z_0 is the origin of the local coordinate system
        % the scalar equation of this plane is n_x*(x-x_0)+n_y*(y-y_0)+n_z*(z-z_0)=0
        % when the normal is not perpendicular to the yz plane,
        % let an arbitrary point on this plane p_1 be
        % p_1 = \hat{x}x_1+\hat{y}y_1+\hat{z}z_1
        %     = \hat{x}x_1+\hat{y}(1.1*y0)+\hat{z}(1.1*z0)
        % then x_1 = x_0-(0.1*n_y*y_0+0.1*n_z*z_0)/n_x
        if nhat(1)~=0
            p1 = [p0(1)-(0.1*nhat(2)*p0(2)+0.1*nhat(3)*p0(3))/nhat(1),1.1*p0(2),1.1*p0(3)];
        elseif nhat(2)~=0
            p1 = [1.1*p0(1),p0(2)-(0.1*nhat(1)*p0(1)+0.1*nhat(3)*p0(3))/nhat(2),1.1*p0(3)];
        elseif nhat(3)~=0
            p1 = [1.1*p0(1),1.1*p0(2),p0(3)-(0.1*nhat(1)*p0(1)+0.1*nhat(2)*p0(2))/nhat(3)];
        else
            error('Invalid normal vector nhat = [0,0,0].\n');
        end
        % let's define the plane local coordinate system
        zhat_prime = nhat/norm(nhat); % let the normal be the local z coordinate
        xhat_prime = p1-p0; % the vector pointing to the point we found be \hat{x}^{\prime}
        xhat_prime = xhat_prime/norm(xhat_prime); % we normalize because we are civilized
        yhat_prime = cross(zhat_prime,xhat_prime); % \hat{y}^{\prime} is orthonormal to both \hat{x}^{\prime} and \hat{y}^{\prime}
        yhat_prime = yhat_prime/norm(yhat_prime);
        Ttilde = [xhat_prime;yhat_prime;zhat_prime]; % linear map to rotate from local to global
        %% transform points coordinates from local to global
        points_global = points_local*Ttilde+p0;
   
    end

end

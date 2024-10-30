function [T,Ttilde] = PlaneCoordinatesRotationMatrices2(nhat,vhat)
%% PlaneCoordinatesRotationMatrices()
% the scalar equation of this plane is n_x*(x-x_0)+n_y*(y-y_0)+n_z*(z-z_0)=0
if abs(dot(nhat,vhat))>1e-12 % we need nhat and vhat to be orthogonal
    error('The vectors need to be orthogonal.\n');
end
% let's define the plane local coordinate system
zhat_prime = nhat/norm(nhat); % let the normal be the local z coordinate
xhat_prime = vhat; % the vector pointing to the point we found be \hat{x}^{\prime}
xhat_prime = xhat_prime/norm(xhat_prime); % we normalize because we are civilized
yhat_prime = cross(zhat_prime,xhat_prime); % \hat{y}^{\prime} is orthonormal to both \hat{x}^{\prime} and \hat{y}^{\prime}
yhat_prime = yhat_prime/norm(yhat_prime);
Ttilde = [xhat_prime;yhat_prime;zhat_prime]; % linear map to rotate from local to global for plane 1
T = inv(Ttilde); % linear map to rotate from global to local for plane 1
end % end PlaneCoordinatesRotationMatrices()
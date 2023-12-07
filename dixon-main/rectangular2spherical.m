function [Spherical_x, Spherical_y, Spherical_z] = rectangular2spherical(theta, phi, z)
%因为画图时使用surf来画图，所以要保证theta,phi,z的矩阵大小相同
%相当于三维空间中，一个点的三个坐标。
R = abs(z);
%z = z - max(z); %如果z小于零，可以只取较大的数
%R=z(z + 20 > 0); 
Spherical_x = R .* sin(theta) .* cos(phi);
Spherical_y = R .* sin(theta) .* sin(phi);
Spherical_z = R .* cos(theta);
end

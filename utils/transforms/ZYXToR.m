function rot = ZYXToR(angles)
% rot = ZYXToR(angles)
% convert ZYX Euler angles to rotation matrix

if numel(angles) ~= 3
    error('numel(angles) must equal 3');
end

phi = angles(1);
theta = angles(2);
psi = angles(3);

rot(1, 1) = cos(theta)*cos(psi);
rot(1, 2) = cos(psi)*sin(theta)*sin(phi) - cos(phi)*sin(psi);
rot(1, 3) = cos(phi)*cos(psi)*sin(theta) + sin(phi)*sin(psi);

rot(2, 1) = cos(theta)*sin(psi);
rot(2, 2) = cos(phi)*cos(psi) + sin(theta)*sin(phi)*sin(psi);
rot(2, 3) = -cos(psi)*sin(phi) + cos(phi)*sin(theta)*sin(psi);

rot(3, 1) = -sin(theta);
rot(3, 2) = cos(theta)*sin(phi);
rot(3, 3) = cos(theta)*cos(phi);

end

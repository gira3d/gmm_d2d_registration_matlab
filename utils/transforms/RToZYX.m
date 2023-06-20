function anglevec = RToZYX(rot)
% anglevec = RToZYX(rot)
% convert rotation matrix to ZYX Euler angles

theta = -asin(rot(3,1));
phi = 0;
psi = 0;

if cos(theta) > eps
    phi = atan2(rot(3,2),rot(3,3));
    psi = atan2(rot(2,1),rot(1,1));
end

anglevec = [phi theta psi]';

end

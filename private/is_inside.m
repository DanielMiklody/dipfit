function inside= is_inside(pos,bnd)
%IS_INSIDE Summary of this function goes here
%   Detailed explanation goes here
inside = bounding_mesh(pos,  bnd.pos, bnd.tri)==1;
end


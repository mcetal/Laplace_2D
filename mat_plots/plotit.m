close all
clear

xgrid;
xg = a;
ygrid;
yg = a;
igrid;
ig = a;
ugrid;
ug = a;

figure(1)
mesh(xg,yg,ig)


figure(2)
geo_contours
target_points

figure(3)
du = ulim(2) - ulim(1);
nc = 40;
vc = ulim(1):du/nc:ulim(2);
contour(xg, yg, ug, vc)
hold on
geo_contours
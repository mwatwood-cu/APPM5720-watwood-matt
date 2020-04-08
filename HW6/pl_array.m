function void = pl_array(px,py)

a = dir('u0*.txt');
ax = dir('x0*.txt');
ay = dir('y0*.txt');


ii = 1;
U = [];
X = [];
Y = [];
for j = 1:py
    uu = []; xx = []; yy = [];
    for i = 1:px
        u = load(a(ii).name);
        x = load(ax(ii).name);
        y = load(ay(ii).name);
        ii = ii + 1;
        uu = [uu u]; xx = [xx x]; yy = [yy y];
    end
    U = [U;uu]; X = [X;xx]; Y = [Y;yy];
end

surf(X,Y,U) 
shading interp
view(30,30)
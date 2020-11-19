function u=fill_lower_triangle(u)
[~,~,nx,ny,nz]=size(u);
for i=1:nx
    for j=1:ny
        for k=1:nz
            for I=1:6
                for J=I:6
                    u(J,I,i,j,k)=u(I,J,i,j,k);
                end
            end
        end
    end
end
end
function [A,Dt2n,R]=form_full_matrix_trimesh_with_new_diffusion_matrix...
    (p,t,HpPhase_nodes,HpPhaseDerx_elements,HpPhaseDery_elements,HpPhaseLaplacian_elements,cdiff,conv_ceff)
%This module uses the elementwise values of derivatives and Laplacian
nt=length(t);np=length(p);
% ar=pdetrg(p,t);%areas of triangles

% Find coefficient matrix
COEF=zeros(3,3,nt);
for k=1:nt
    K = [1,p(1,t(1,k)),p(2,t(1,k));1,p(1,t(2,k)),p(2,t(2,k));1,p(1,t(3,k)),p(2,t(3,k))];
    COEF(:,:,k) = inv(K);
end
% figure(1); imshow(abs(COEF(:,:,1)),[]);title('COEF')

% Find Convection matrix
C=sparse(nt,np);
for k=1:nt
    C(k,t(1:3,k)) = COEF(2,:,k)*HpPhaseDerx_elements(k)+COEF(3,:,k)*HpPhaseDery_elements(k);
end
% figure(2); imshow(abs(full(C)),[]);title('C')

% Find Reaction matrix
R=sparse(nt,np);
for k=1:nt
    R(k,t(1:3,k)) = HpPhaseLaplacian_elements(k)/3;
end
% figure(3); imshow(abs(full(R)),[]);title('R')

%Find Diffusion matrix
D=sparse(nt,np);
%gradient matrices from nodes to triangles

Cx=sparse(nt,np);Cy=sparse(nt,np);
for k=1:nt
    Cx(k,t(1:3,k)) = COEF(2,:,k);
    Cy(k,t(1:3,k)) = COEF(3,:,k);
end
% figure(4); imshow(abs(full(Cx)),[]);title('Cx')
% figure(5); imshow(abs(full(Cy)),[]);title('Cy')

    %interpolation matrix from triangles to nodes with triangle areas not considered
    Dt2n=sparse(np,nt);
    for i=1:np
        [rows,cols,values]=find(t(1:3,:)==i);%alternatively use pdeneigh
        for k=1:length(cols)  
            Dt2n(i,cols(k))=1/length(cols);
        end
    end    
 
%     figure(6); imshow(abs(full(Dt2n)),[]);title('Dt2n')


D=Cx*Dt2n*Cx+Cy*Dt2n*Cy;
% figure(7); imshow(abs(full(D)),[]);title('D')

A=(conv_ceff*C+R+cdiff*D);


%some usefull codes
% %interpolation matrix from triangles to nodes with triangle areas considered
% Dt2n=sparse(np,nt);
% for i=1:np
%     ii=find(ts(1:3,:)==i);iii=floor((ii-1)/3)+1;iiii=unique(iii);%alternatvely use pdeneigh
%     for k=1:length(iiii)  Dt2n(i,iiii(k))=ar(iiii(k))/sum(ar(iiii));end
% end

% for i=1:np
%     ii=find(ts(1:3,:)==i);iii=floor((ii-1)/3)+1;iiii=unique(iii);%alternatvely use pdeneigh
%     for k=1:length(iiii)  Dt2n(i,iiii(k))=1/length(iiii);end
% end
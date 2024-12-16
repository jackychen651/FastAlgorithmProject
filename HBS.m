clear;clc;close all;
%h = 2.0^(-5);
omega = 25.0;
L=8;
N=2^L;
h=1/N;
a = 1;
b = 1;

x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

[X,Y] = meshgrid(x,y);

X_vec = X(:);
Y_vec = Y(:);
tic
G=zeros(N^2, N^2);
% for i=1:N^2
%     for j=1:N^2
%         if j ~= i
%             G(i,j) = h^2 *1j / 4* besselh(0, 1, omega * sqrt((X_vec(i)-X_vec(j))^2+(Y_vec(i)-Y_vec(j))^2));
%         end
%     end
% end

X_ = repmat(X_vec, 1, N^2);
Y_ = repmat(Y_vec, 1, N^2);

dists = sqrt((X_ - X_.') .^ 2 + (Y_ - Y_.') .^ 2);
dists(1:N^2+1:end) = 1;
B_ = besselh(0, 1, omega * dists);
B_(1:N^2+1:end) = 0;

G = (h^2 * 1j / 4) * B_;

nu = @(x,y) -1.5*exp(-160*(x.^2 + y.^2));

B=omega^2*nu(X_vec,Y_vec);
B = diag(B);
f = -omega^2*nu(X_vec,Y_vec).*exp(omega*1i*X_vec);
toc

disp('All matrices formed successfully');

tic
sigma = (eye(N^2)+B*G)\f;


u = G*sigma;
u_inc = exp(omega*1i*X);
u = reshape(u, N, N);
toc

M = nu(X,Y);

figure(1); clf();
imagesc(M);
colorbar;
filename = sprintf('Figure1_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);

figure(2); clf();
imagesc(real(reshape(u, n, m)));
colorbar;
filename = sprintf('Figure2_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);

figure(3); clf();
imagesc(real(reshape(u + u_inc, n, m)));
colorbar;
filename = sprintf('Figure3_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);


%%


%first find U and V
depth=L;
Index=cell(depth,2^(depth-1));
for l=1:depth
    for j=1:2^(l-1)
        Index{l,j}=(1+(j-1)*2^(2*L-l+1):j*2^(2*L-l+1));
    end
end

U = 1:(2^(2*L));
Index_c = cell(depth, 2^(depth-1));
for l = 1:depth
    for j = 1:2^(l-1)
        Index_c{l, j} = setdiff(U, Index{l, j});
    end
end


k=4;
%[J_row,U]=id_row(G(Index{L-1,1},Index{L-1,2}),2);
U=cell(depth,2^(depth-1));
V=cell(depth,2^(depth-1));
W=cell(depth,2^(depth-1));
J=cell(depth,2^(depth-1));
skel_id=cell(depth,2^(depth-1));
J_row=cell(depth,2^(depth-1));
J_col=cell(depth,2^(depth-1));
skel_row_id=cell(depth,2^(depth-1));
skel_col_id=cell(depth,2^(depth-1));
for l=2:depth
    for j=1:2^(l-1)
        [U{l,j},J_row{l,j}]=id_row(G(Index{l,j}, Index_c{l, j}),k);
        [V{l,j},J_col{l,j}]=id_row(G(Index_c{l, j}, Index{l,j}).',k);
        skel_row_id{l,j}=Index{l,j}(J_row{l,j});
        skel_col_id{l,j}=Index{l,j}(J_col{l,j});
        W{l,j}=[G(Index{l,j}, Index_c{l, j}), G(Index_c{l, j}, Index{l,j}).'];
        [U{l,j},J{l,j}]=id_row(W{l,j},k);
        skel_id{l,j}=Index{l,j}(J{l,j});
    end
end
X=cell(depth,2^(depth-1));
S=cell(depth,2^(depth-1));
for l=depth:depth
    for j=1:2^(l-1)
        B_tmp=B(Index{l,j}, Index{l,j});
        G_tmp=G(Index{l,j},Index{l,j});
        X{l,j}=inv(eye(size(Index{l,j}, 2))+B_tmp*G_tmp);
        S{l,j}=V{l,j}.'*X{l,j}*B_tmp*U{l,j};
    end
end
for l=depth-1:-1:1
    for j=1:2^(l-1)
        alpha=Index{l+1,2*j-1};
        beta=Index{l+1,2*j};
        X{l,j}=inv([eye(size(alpha, 2)),S{l+1,2*j-1}*G(alpha, beta);S(l+1,2*j)*G(beta, alpha), eye(size(alpha, 1))]);
        if l~=1
            S{l,j}=V(l,j).'*X(l,j)*[S(l+1,2*j-1),zeros(size(alpha, 1));zeros(size(alpha, 1)),S(l+1,2*j)]*U(l,j);
        end
    end
end


function [X,I_s]=id_row(A,k)
    m=size(A,1);
    [~,S,J]=qr(A',0);
    T=(S(1:k,1:k))\S(1:k,(k+1):m);
    X=zeros(m,k);
    X(J,:)=[eye(k),T]';
    I_s=J(1:k);
end
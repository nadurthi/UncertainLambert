function pert_transorbit=convert_transferorbit_mat2struct_perturbed(M)


% M=[pert_transorbit.N,pert_transorbit.delV,pert_transorbit.delVD,pert_transorbit.delVA,pert_transorbit.r1(:)',pert_transorbit.r2(:)',pert_transorbit.v1(:)',pert_transorbit.v2(:)',...
%     pert_transorbit.v0A(:)',pert_transorbit.v0B(:)'];
% M=[M,pert_transorbit.true_timevec(:)',pert_transorbit.true_Altitude(:)',reshape(pert_transorbit.true_traj,1,size(pert_transorbit.true_traj,1)*6)];


pert_transorbit.N=M(1);
pert_transorbit.delV=M(2);
pert_transorbit.delVD=M(3);
pert_transorbit.delVA=M(4);
pert_transorbit.r1=M([5,6,7]);
pert_transorbit.r2=M([8,9,10]);
pert_transorbit.v1=M([11,12,13]);
pert_transorbit.v2=M([14,15,16]);
pert_transorbit.v0A=M([17,18,19]);
pert_transorbit.v0B=M([20,21,22]);
m=length(M(23:end));
n=m/8;
pert_transorbit.timevec=M(23:23+n-1);
pert_transorbit.true_Altitude=M(23+n:23+n+n-1);
pert_transorbit.traj=reshape(M(23+2*n:end),n,6 );

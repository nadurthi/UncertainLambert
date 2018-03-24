function M=convert_transferorbit_struct2mat(transorbit)
% as cells and structs use a lot of memory, they will be converted to a
% vector. a corresposnding function convert_transferorbit_mat2struct can be
% used to get the struct back
% - convert_transferorbit_struct2mat and convert_transferorbit_mat2struct go
%   hand in hand

if strcmp(transorbit.branch,'lower')
    transorbit.branch=-1;
else
    transorbit.branch=1;
end

M=[transorbit.N,transorbit.branch,transorbit.a,transorbit.s,transorbit.c,transorbit.theta,transorbit.t_des,transorbit.a_m,transorbit.delV,transorbit.delVD,transorbit.delVA,transorbit.r1(:)',transorbit.r2(:)',transorbit.v1(:)',transorbit.v2(:)',transorbit.v0A(:)',transorbit.v0B(:)'];
M=[M,transorbit.normalized_transfer_timesteps(:)',transorbit.true_Altitude(:)',reshape(transorbit.traj,1,size(transorbit.traj,1)*6)];
    
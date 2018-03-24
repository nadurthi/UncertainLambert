function transorbit=convert_transferorbit_mat2struct(M)
% - convert M , a vector, into a struct with all the required parameters of a
%   transferorbit structure
% - convert_transferorbit_struct2mat and convert_transferorbit_mat2struct go
%   hand in hand


transorbit.N=M(1);
if M(2)==-1
    transorbit.branch='lower';
else
    transorbit.branch='upper';
end
transorbit.a=M(3);
transorbit.s=M(4);
transorbit.c=M(5);
transorbit.theta=M(6);
transorbit.t_des=M(7);
transorbit.a_m=M(8);
transorbit.delV=M(9);
transorbit.delVD=M(10);
transorbit.delVA=M(11);
transorbit.r1=M([12,13,14]);
transorbit.r2=M([15,16,17]);
transorbit.v1=M([18,19,20]);
transorbit.v2=M([21,22,23]);
transorbit.v0A=M([24,25,26]);
transorbit.v0B=M([27,28,29]);
m=length(M(30:end));
n=m/8;
transorbit.normalized_transfer_timesteps=M(30:30+n-1);
transorbit.true_Altitude=M(30+n:30+n+n-1);
transorbit.traj=reshape(M(30+2*n:end),n,6 );


    
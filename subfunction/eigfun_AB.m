function [finteVf,finiteLam,A,B]=eigfun_AB(r,D,N,w,m,Ratio,Mx,M_theta,rou0,P0,c02,Boundary,z_t,z_h)
    gama=1.4;Dr_rou0=D*rou0;Dr_Mx=D*Mx;Dr_M_theta=D*M_theta;
    
%    set up matrices A and B  
    A21=[];
    for ii=1:N+1
        A1{ii}=[-sqrt(-1)*w+M_theta(ii)*sqrt(-1)*m/r(ii)  rou0(ii)*M_theta(ii)^2/r(ii)+gama*P0(ii)/r(ii)  gama*P0(ii)*sqrt(-1)*m/r(ii)              0                                0;...
                    0                                     -sqrt(-1)*w+M_theta(ii)*sqrt(-1)*m/r(ii)         -2*M_theta(ii)/r(ii)                      0                -M_theta(ii)^2/rou0(ii)/r(ii);...
                 sqrt(-1)*m/r(ii)/rou0(ii)                  M_theta(ii)/r(ii)+Dr_M_theta(ii)               -sqrt(-1)*w+M_theta(ii)*sqrt(-1)*m/r(ii)   0                              0;...
                    0                                         Dr_Mx(ii)                                         0                        -sqrt(-1)*w+M_theta(ii)*sqrt(-1)*m/r(ii)    0;...
                    0                                         Dr_rou0(ii)+rou0(ii)/r(ii)                       rou0(ii)*sqrt(-1)*m/r(ii)              0               -sqrt(-1)*w+M_theta(ii)*sqrt(-1)*m/r(ii)];
        A21=[A21;gama*P0(ii);1/rou0(ii);0;0;rou0(ii)];
    
        B{ii}=[Mx(ii)  0     0  gama*(P0(ii)) 0;
                0     Mx(ii)  0      0      0;
                 0       0   Mx(ii)  0      0;
             1/rou0(ii)  0     0   Mx(ii)   0;
                0        0     0   rou0(ii) Mx(ii);];
    end    
    A1 = blkdiag(A1{:});B = blkdiag(B{:});
    A2=A21.*(kron(D,full(sparse(1,2,1,5,5)))+kron(D,full(sparse(2,1,1,5,5)))+kron(D,full(sparse(5,2,1,5,5))));

% set up boundary for matrices A and B
    
    switch Boundary
        case {1} %Hard inside and hard outside
            A2([1 5],2)=0;%!!!
            A1([2],1:5)=[0 1 0 0 0];   
            A1([1:5],2)=[0;1;0;0;0]; 
            B(2,2)=0;
            A2([end-4 end],end-3)=0;%!!!
            A1([end-3],end-4:end)=[0 1 0 0 0];   
            A1([end-4:end],end-3)=[0;1;0;0;0]; 
            B(end-3,end-3)=0;           
        case{2} %Hard inside and soft outside
            A1(4,1:5)=[-sqrt(-1)*w+M_theta(1)*sqrt(-1)*m/r(1) sqrt(-1)*z_t*w 0 0 0;];
            B(4,1:5)=[Mx(1) 0 0 0 0];  
            A1([end-2],end-3:end)=[0 1 0 0;];A2([end-3],2)=0;%!!  
            A1(end,end-2)=0;A1(end-1,end-2)=0;A1(end-3,end-2)=0;             
            B(end-2,end-2)=0;
        case{3} %Soft inside and hard outside
            A1(end-1,end-4:end)=[-sqrt(-1)*w+M_theta(end)*sqrt(-1)*m/r(end) -sqrt(-1)*z_t*w 0 0 0;];
            B(end-1,end-4:end)=[Mx(end) 0 0 0 0];   
            A2([1 5],2)=0;%!!!
            A1([2],1:5)=[0 1 0 0 0];   
            A1([1:5],2)=[0;1;0;0;0]; 
            B(2,2)=0;  
        case{4} %Soft inside and soft outside
            A1(4,1:5)=[-sqrt(-1)*w+M_theta(1)*sqrt(-1)*m/r(1) sqrt(-1)*z_t*w 0 0 0;];
            B(4,1:5)=[Mx(1) 0 0 0 0];    
            A1(end-1,end-4:end)=[-sqrt(-1)*w+M_theta(end)*sqrt(-1)*m/r(end) -sqrt(-1)*z_t*w 0 0 0;];
            B(end-1,end-4:end)=[Mx(end) 0 0 0 0];  
    end
        
        A=A1+A2;
        [V,Lam]=eig((-A-B),(B-A));
        lam = sqrt(-1)*(1+2./(diag(Lam)-1))/r(1);

        for kk=1:5*(N+1)
            vf{kk} = [reshape(V(:,kk),[5,N+1])];
        end
        finiteMode= find(real(lam)>-70000&real(lam)<70000&imag(lam)>-70000&imag(lam)<70000&(real(lam).^2+imag(lam).^2)>1E-3);
        finiteLam=lam(finiteMode);
        for kk=1:length(finiteMode)
            finteVf{kk}=vf{finiteMode(kk)};
        end

end
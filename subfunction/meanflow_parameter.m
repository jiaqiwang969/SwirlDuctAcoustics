 function  [rou0_3D,P0_3D,s0_3D,Mx_3D,M_theta_3D]=meanflow_parameter(rou0,P0,s0,Mx,M_theta,x_pole,thetaNumber)

        kkk=1;
        for xk0=x_pole
            rou0_3D(:,:,kkk)=rou0*ones(1,thetaNumber);
            P0_3D(:,:,kkk)=P0*ones(1,thetaNumber);
            rou0_3D(:,:,kkk)=rou0*ones(1,thetaNumber);
            s0_3D(:,:,kkk)=s0*ones(1,thetaNumber);
            Mx_3D(:,:,kkk)=Mx*ones(1,thetaNumber);
            M_theta_3D(:,:,kkk)=M_theta*ones(1,thetaNumber);
            %Mr_3D(:,:,kkk)=zeros(length(rou0),1)*ones(1,360);
            kkk=kkk+1;
        end 


    end




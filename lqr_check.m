function [A,B,cost,u,x_0,T,q,r]=lqr_check
%optimal deterministic LQR control for mass spring damper with monte carlo simulation
clc
clear
close all

global t N_u N_opt r x_0 dt eps K q r1 A B x u cost


T=8; %final time
q=[2 0;0 2];%2*eye(2); %weight for state in stochastic cost function
r=1; %weight for control in stochastic cost function

N_u=4000; %No. of time samples for state propagation
N_opt=100; %No. of time samples for control optimization
r1=N_u/N_opt; %ratio of sampling

t=linspace(0,T,N_u+1); %simulation time vector
dt=t(2)-t(1);
eps=.1;


%system matrices
m=1;
k=1;
c=2;
A=[0 1;-k/m -c/m];
B=[0 1]';
[KK,S,e]=lqr(A,B,q,r,zeros(2,1)); %u=-KK x// S: solution of ricatti equation//e=eig(A-B KK)
% keyboard
x=zeros(2,N_u+1); %N_u+1 by 1
x_0=[1 0]'; %initial state
x(:,1)=x_0;
xdot=zeros(2,N_u);
J(1)=0;

u=zeros(1,N_u+1);
for kk=1:N_u
   u(1,kk)=-KK*x(:,kk);
   x(:,kk+1)=x(:,kk)+dt*(A*x(:,kk)+B*u(1,kk));
   xdot(:,kk)=A*x(:,kk)+B*u(1,kk);
   integrand=q(1,1)*x(1,1:kk).^2+q(2,2)*x(2,1:kk).^2+r*u(1,1:kk).^2;
   if kk>=2
       J(kk)=trapz(t(1:kk),integrand);
   end
end
u(1,N_u+1)=-KK*x(:,N_u+1);

figure(1);plot(t,x(1,:),t,x(2,:),t,u,'--');
legend('$x^{*}_1$ optimal trajectory','$x^{*}_2$ optimal trajectory','$u^{*}$ optimal control','Interpreter','Latex');
grid on;

figure(2);plot(t(2:end),J);
legend('J: cost function','Interpreter','Latex');
grid on;

cost=J(end);

M=1e6;
K=.15; %channel width
sigma=.05;
% monte_carlo_LQR(M,K,sigma)

end


function monte_carlo_LQR(M,K,sigma)

    global x_0 N_u A B q r dt t x u cost
%     sample_accu_cost=zeros(N_u+1,M);
    %derive M monte carlo samples 
    x_eps1=zeros(2,N_u+1);
    x_eps1(:,1)=x_0;
    sum_path=zeros(2,N_u+1);
    sample_cost=zeros(M,1);
    first_exit=zeros(M,1);
    temp_sum=zeros(2,N_u+1);
    temp_sum2=0;
    temp_sum3=zeros(2,N_u+1);
    counter=0;
    
    sample_accu_cost=zeros(1,N_u+1);
    save_traj=[];
    save_cost=[];
    save_error=[];

    for m=1:M

        dW=sigma*sqrt(dt)*randn(2,N_u+1);
        for k=1:(N_u)
            x_eps1(:,k+1)=x_eps1(:,k)+dt*(A*x_eps1(:,k)+B*u(k))+dW(:,k);
            sample_accu_cost(k+1)=sample_accu_cost(k)...
                +dt*(q(1,1)*x_eps1(1,k+1)^2+q(2,2)*x_eps1(2,k+1)^2+r*u(k)^2);
        end
        sample_cost(m)=trapz(t,q(1,1)*x_eps1(1,:).^2+q(2,2)*x_eps1(2,:).^2+r*u.^2);

        err=x_eps1-x; %e_epsilon in notation table 7.3
        temp_sum=temp_sum+err;
        
        if rem(m,10000) == 0
           disp(['MonteCarlo for LQR sim. # ',num2str(m)]) 
           save_traj=[save_traj; x_eps1];
           save_error=[save_error; err];
           save_cost=[save_cost;sample_accu_cost];
        end
        
        

        temp_sum2=temp_sum2+max(max(err));
        temp_sum3=temp_sum3+err.^2;
        [row,col]=find(abs(err)>=K,1,'first');
        if col
            first_exit(m)=(col-1)*dt;
            counter=counter+1;
        end
        
        sum_path=sum_path+x_eps1;

    end

    mean_path=sum_path/M;
    var_error=temp_sum3/M;
    
    E_e_eps=temp_sum/M;%mean(e_eps); %expected(e_epsilon) average error
    int_E_e_eps1=trapz(t,E_e_eps(1,:)); %integral {expected(e_epsilon1)}dt
    int_E_e_eps2=trapz(t,E_e_eps(2,:)); %integral {expected(e_epsilon2)}dt

    max_E_e_eps=max(max(E_e_eps)); %max{expected(e_epsilon)}
    [row,col]=find(E_e_eps==max_E_e_eps); %t{ max{expected(e_epsilon)} }
    tmax_e_eps=t(col);

    E_max_eeps=temp_sum2/M;  %expected[ max(e_epsilon) ]

    exit_probability=counter/M %based on trajectory exiting the channel around the nominal value

    E_D=mean(sample_cost); %expexted(D) or expected(v) table 7.3 page 261
    var_D=var(sample_cost); %variance(D)
    
    err_cost=sample_cost-cost;
    col=find(abs(err_cost)>=K);
    exit_prob_cost=length(col)/M;
    
    save('lqr_mc.mat','exit_probability','E_max_eeps',...
        'max_E_e_eps','tmax_e_eps','E_e_eps','int_E_e_eps1',...
        'int_E_e_eps2','E_D','var_D','exit_prob_cost'...
        ,'mean_path','var_error','save_cost','save_traj');
    
end

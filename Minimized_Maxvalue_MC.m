%该程序包含了一个基于蒙特卡洛法的，在多频合成时，通过设置各频率分量相位使得合成波表的最大值最小(峰均比最小)的算法
clear

Sa=800;%采样点数
fs=200;%采样率
t=(0:Sa-1)/fs;%time
%f_list=[66e6 67e6 68e6 69e6 70e6];%频谱量
%f_list=(70:1:79)*1e6;
f_list=10:19;
fN=length(f_list);
%a_list=[0.9 0.95 1.15 1.15 0.95];%频率对应幅值
a_list=ones(1,fN);
sum_a=sum(a_list);
N=2e3;%蒙特卡洛迭代次数
n=N/10;%每次查漏迭代次数，少一个量级
temp=1000;

%蒙特卡洛搜寻核心
%第一次粗略搜寻，在整个2pi域内搜索出一组较好的种子相位
for i=1:N*2
    phi_list=2*pi*rand(1,fN-1);
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%在第一组种子相位上减少搜索范围到pi，进行第一次精细搜索
old_phi=temp_phi;
for i=1:N
    phi_list=old_phi+pi*rand(1,fN-1)-pi/2;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%在上一组最优相位基础上，减少搜索范围到pi/4，进行第二次精细搜索
old_phi=temp_phi;
for i=1:N
    phi_list=old_phi+pi/4*rand(1,fN-1)-pi/8;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%避免边界还有更优的解，在上一组最优相位基础上，保持范围pi/4，查漏
old_phi=temp_phi;
for i=1:n
    phi_list=old_phi+pi/4*rand(1,fN-1)-pi/8;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%在上一组最优相位基础上，减少搜索范围到pi/16，进行第三次精细搜索
old_phi=temp_phi;
for i=1:N
    phi_list=old_phi+pi/16*rand(1,fN-1)-pi/32;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%避免边界还有更优的解，在上一组最优相位基础上，保持范围pi/16，查漏
old_phi=temp_phi;
for i=1:n
    phi_list=old_phi+pi/16*rand(1,fN-1)-pi/32;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%在上一组最优相位基础上，减少搜索范围到pi/64，进行第四次精细搜索
old_phi=temp_phi;
for i=1:N
    phi_list=old_phi+pi/64*rand(1,fN-1)-pi/128;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end
%避免边界还有更优的解，在上一组最优相位基础上，保持范围pi/64，查漏
old_phi=temp_phi;
for i=1:n
    phi_list=old_phi+pi/64*rand(1,fN-1)-pi/128;
    sum_s=a_list(1).*sin(2*pi*f_list(1)*t);
    for j=1:fN-1
        sum_s=sum_s+a_list(j+1).*sin(2*pi*f_list(j+1)*t+phi_list(j));
    end
    max_s=max(sum_s);
    if max_s<temp
        temp=max_s;
        temp_phi=phi_list;
        temp_sum=sum_s;
    end
end


%未优化叠加
sum_0=a_list(1).*sin(2*pi*f_list(1)*t);
for j=1:fN-1
    sum_0=sum_0+a_list(j+1).*sin(2*pi*f_list(j+1)*t);
end

temp_phi=[0 temp_phi];
for i=1:fN
    if i==1
        str=append(num2str(a_list(i)),"*sin(2*pi*",num2str(f_list(i)),"*x)");
    else
        str=append(str,"+",num2str(a_list(i)),"*sin(2*pi*",num2str(f_list(i)),"*x+",num2str(temp_phi(i),'%.2f'),")");
    end
end
str=append("0.1*(",str,")");
disp(str);
PAR=10*log10(temp.^2/sum(a_list.^2));
disp(append("PAR=",num2str(PAR),"dB"));

figure(1)
subplot(2,1,1)
plot(t,sum_0)
title(["未优化" append("PAR=",num2str(10*log10(max(sum_0).^2/sum(a_list.^2))),"dB")])
xlabel("t/s")
ylabel("A")
ylim([-sum_a sum_a])
subplot(2,1,2)
plot(t,temp_sum)
title(["相位优化后" append("PAR=",num2str(PAR),"dB")])
xlabel("t/s")
ylabel("A")
ylim([-sum_a sum_a])

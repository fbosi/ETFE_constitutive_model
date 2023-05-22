%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads images
%
% Written by: A. Comitti, 2023 ~  a.comitti@ucl.ac.uk
%
% Copyright F. Bosi, 2023; A. Comitti, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
path = pwd;
%%%%%%%%%%%% loading the model constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = pwd;
opt1 = readmatrix(strcat(path,'\D11_Opt_1e5Mpa.csv'));
opt2 = readmatrix(strcat(path,'\D22_Opt_1e5Mpa.csv'));
opt6 = readmatrix(strcat(path,'\D66_Opt_1e5Mpa.csv'));
ea11 = readmatrix(strcat(path,'\EA11_1e5Mpa.csv'));    
ea22 = readmatrix(strcat(path,'\EA22_1e5Mpa.csv'));
ea66 = readmatrix(strcat(path,'\EA66_1e5Mpa.csv'));
ea11 = ea11/3 + ea22/3 + ea66/3;
factor = 2.303*8.31446261815324;
TRef = 20 + 273.15;
poisson = 0.43;
thickness = 0.2;

%%%%%%%%% TEST DATA LOADING; comment out the test not needed %%%%%%%%%%%%%
%%%% the raw data are all expressed in engineering quantities %%%%%%%%%%%%




%%%% uniaxial constant strain rate tests data from Instron, MD %%%%%%%%%%%
% rawData = rmmissing(readmatrix(strcat(path,'\OUTRaw_MD_CSR_01_T40.csv')));
% tempData = rmmissing(readmatrix( ...
%     strcat(path,'\OUTTemp_MD_CSR_01_T40.csv')));
% % temperature vector synchronisation
% for i = 1: size(tempData,2)-1
%     T{i} = rmmissing(interp1(tempData(:,1)-1,tempData(:,i+1), ...
%         rawData(:,1+15*(i-1)),'linear','extrap'))+273.15; 
% end
% 
% % interpolation on a strain vector
% e11 = transpose(0:0.001:0.02);
% for i = 1: size(tempData,2)-1
%     [e11C{i}, index,~] = unique(rawData(:,2+15*(i-1)),'stable');
%     trefC(:,i) = interp1(e11C{i},rawData(index,1+15*(i-1)),e11,'linear');
%     e22C(:,i) =  interp1(e11C{i},rawData(index,5+15*(i-1)),e11,'linear');
%     e12C(:,i) = zeros(length(e11),1);
%     s11C(:,i) = interp1(e11C{i},rawData(index,4+15*(i-1)),e11,'linear');
%     s22C(:,i) = zeros(length(e11),1);
%     s12C(:,i) = zeros(length(e11),1);
%     tempC(:,i) = interp1(e11C{i},T{i}(index),e11,'linear');
% end
% tref = mean(trefC,2,'omitnan'); 
% e22 = mean(e22C,2,'omitnan'); 
% e12 = zeros(length(e11),1);
% s11 = mean(s11C,2,'omitnan'); 
% s22 = zeros(length(e11),1);
% s12 = zeros(length(e11),1);
% temp = mean(tempC,2,'omitnan'); 
% s11Sd =  std(s11C,0,2,'omitnan'); 
% e11Sd =  zeros(length(e11),1);
% E = 1030.02;



%%%% uniaxial relaxation tests data from Instron, MD %%%%%%%%%%%%%%%%%%%%%
rawData = rmmissing(readmatrix(strcat( ...
    path,'\OUTRaw_MD_Rel_9MPa_T20.csv')));
tempData = rmmissing(readmatrix(strcat( ...
    path,'\OUTTemp_MD_Rel_9MPa_T20.csv')));
% temperature vector synchronisation
for i = 1: size(tempData,2)-1
    T{i} = rmmissing(interp1(tempData(:,1),tempData(:,i+1), ...
        rawData(:,1+9*(i-1)),'linear',tempData(end,i+1)))+273.15; 
end

% interpolation on a time vector
tref = transpose(logspace(-1,4,40000));
for i = 1: size(tempData,2)-1
    [t{i},index,~] = unique(rawData(:,1+9*(i-1)),'stable');
    e11C(:,i)  =  interp1(t{i}, ...
        smooth(rawData(:,2+9*(i-1)),10),tref,'nearest');
    e22C(:,i) = interp1(t{i},rawData(:,4+9*(i-1)),tref,'nearest');
    e12C(:,i) = zeros(length(tref),1);
    s11C(:,i) = interp1(t{i}, ...
        smooth(rawData(:,3+9*(i-1)),10),tref,'nearest');
    s22C(:,i) = zeros(length(tref),1);
    e12C(:,i) = zeros(length(tref),1);
    tempC(:,i)  = interp1(t{i},T{i}(index),tref,'nearest');
end    
e11 = mean(e11C,2,'omitnan'); 
e22 = mean(e22C,2,'omitnan'); 
e12 = zeros(length(e11),1);
s11 = mean(s11C,2,'omitnan'); 
s22 = zeros(length(e11),1);
s12 = zeros(length(e11),1);
temp = mean(tempC,2,'omitnan'); 
s11Sd =  std(s11C,0,2,'omitnan');
e11Sd =  std(e11C,0,2,'omitnan');
E = 1219;

%%%%%%%%%%%%% Recursive algorithm integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for q = 1:length(opt1)-1
            tau(q)=opt1(q+1,2);
            D110=opt1(1,1);
            D11j(q)=opt1(q+1,1);
            D220=opt2(1,1); 
            D22j(q)=opt2(q+1,1);
            D120=-poisson*D110;
            D12j(q)=-poisson*D11j(q);
            D210= D120; 
            D21j(q)=D12j(q);
            D130 =  D120;
            D13j(q) = D12j(q);
            D230 =  D120;
            D23j(q) = D12j(q);
            D660 =  opt6(1,1);
            D66j(q) = opt6(q+1,1);
end

%%%%% initialisation
s11_R = zeros(length(e11),1) ;% engineering stress 11  
s22_R = zeros(length(e11),1) ;% engineering stress 22
s12_R = zeros(length(e11),1) ;% engineering stress 12
e33_R = zeros(length(e11),1); % out of plane engineering strain
q11 = zeros(length(tau),1) ;
q11old = zeros(length(tau),1) ;
q12 = zeros(length(tau),1) ;
q12old = zeros(length(tau),1) ;
q21 = zeros(length(tau),1) ;
q21old = zeros(length(tau),1) ;
q22 = zeros(length(tau),1) ;
q22old = zeros(length(tau),1) ;
q13 = zeros(length(tau),1) ;
q13old = zeros(length(tau),1) ;
q23 = zeros(length(tau),1) ;
q23old = zeros(length(tau),1) ;
q66 = zeros(length(tau),1) ;
q66old = zeros(length(tau),1) ;

%starting the thing
for n = 2:1:length(e11)
        %Shift factor and reduced time
        aT = 10^(-ea11/factor*(-1/temp(n) + 1/TRef));
        dtr = (tref(n)-tref(n-1))/(aT);
        %Initialization of parameters
        F11 = 0;        F11a = 0;
        D11 = D110;
        F12 = 0;        F12a = 0;
        D12 = D120;
        F21 = 0;        F21a = 0;
        D21 = D210;
        F22 = 0;        F22a = 0;
        D22 = D220;
        F13 = 0;
        D13 = D130;
        F23 = 0;
        D23 = D230; 
        F66 = 0;
        D66 = D660; 
        %hereditary storage 
        q11old = q11;
        q12old = q12;
        q21old = q21;
        q22old = q22;
        q13old = q13;
        q23old = q23;
        q66old = q66;
        %looping through the Prony series
        for c = 1:length(tau)
            dtrexp = dtr/tau(c);
            if dtrexp<1e-10
                g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
            else
                g = (1-exp(-dtrexp))/dtrexp;
            end
            D11 = D11 + D11j(c)*(1-g);
            D12 = D12 + D12j(c)*(1-g);
            D22 = D22 + D22j(c)*(1-g);
            D21 = D21 + D21j(c)*(1-g);
            D13 = D13 + D13j(c)*(1-g);
            D23 = D23 + D23j(c)*(1-g);
            D66 = D66+ D66j(c)*(1-g);
            F11 = F11 + D11j(c)*(exp(-dtrexp)*q11old(c)-g*s11_R(n-1));
            F12 = F12 + D12j(c)*(exp(-dtrexp)*q12old(c)-g*s11_R(n-1));
            F21 = F21 + D21j(c)*(exp(-dtrexp)*q21old(c)-g*s22_R(n-1));
            F22 = F22 + D22j(c)*(exp(-dtrexp)*q22old(c)-g*s22_R(n-1));
            F13 = F13 + D13j(c)*(exp(-dtrexp)*q13old(c)-g*s11_R(n-1));
            F23 = F23 + D23j(c)*(exp(-dtrexp)*q23old(c)-g*s22_R(n-1));
            F66 = F66 + D66j(c)*(exp(-dtrexp)*q66old(c)-g*s12_R(n-1)); 
        end

%calculation of stresses 
        detD = D11*D22-D21*D12;
        s11_R(n) = 1/detD *(D22*(e11(n)+F11+F12)-D21*e22(n)-D12*(F12+F22));
        s22_R(n) = 1/detD *(D11*(e22(n)+F22+F21)-D12*e11(n)-D21*(F21+F11));
        s12_R(n) = (e12(n)+F66)/D66;
        e33_R(n) = D13*s11_R(n) +  D23*s22_R(n); 
        
%calculation of the hereditary integral at the end of the step
        for c = 1:length(tau)
            dtrexp = dtr/tau(c);
            if dtrexp<1e-10
                g = 1-dtrexp/2+(dtrexp)^2/6-(dtrexp)^3/24;
            else
                g = (1-exp(-dtrexp))/dtrexp;
            end
            q11(c) = exp(-dtrexp)*q11old(c)+(s11_R(n)-s11_R(n-1))*g;
            q12(c) = exp(-dtrexp)*q12old(c)+(s11_R(n)-s11_R(n-1))*g;        
            q21(c) = exp(-dtrexp)*q21old(c)+(s22_R(n)-s22_R(n-1))*g;
            q22(c) = exp(-dtrexp)*q22old(c)+(s22_R(n)-s22_R(n-1))*g;
            q13(c) = exp(-dtrexp)*q13old(c)+(s11_R(n)-s11_R(n-1))*g;        
            q23(c) = exp(-dtrexp)*q23old(c)+(s22_R(n)-s22_R(n-1))*g;   
            q66(c) = exp(-dtrexp)*q66old(c)+(s12_R(n)-s12_R(n-1))*g; 
        end
end

%% plot 
fig = figure();  
Col = ['k','r','b','m','g'];
Blue = 1/255*[150,150,245];
Orange = 1/255*[245,150,150];

%%%%% Stress - strain plot %%%%%%
plot(e11,s11,'LineStyle', '-','Color',Blue,'LineWidth',0.5, ...
    'DisplayName', strcat('Experimental'))
hold on
erev = [e11' flip(e11')];
inbetween = [s11-s11Sd;flip(s11Sd+s11)];
h1=  fill(erev, inbetween, Blue,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
    "DisplayName", "Experimental");
h2 =   plot(e11,s11_R,'LineStyle','-.','Color',Col(4),'LineWidth',2, ...
    'DisplayName',strcat('Model'));
h3 =   plot(e11,E*e11,'LineStyle',':','Color',Col(3),'LineWidth',2, ...
    'DisplayName',strcat('Experimental modulus'));
hold off
xlim([0 max(e11)*1.2])
ylim([0 max(s11)*1.2])
set(gca,'fontname','Times New Roman') 
set(gcf, 'color', [1 1 1])
fontsize(gca,16,"points")
fig.Units = 'centimeters';
set(gcf,'Position',[3 3 19 15.2])
xlabel("$\epsilon$",'FontName','Times New Roman','FontSize',18, ...
    'Interpreter','latex')
ylabel("$\sigma$ [MPa]",'FontName','Times New Roman','FontSize',18, ...
    'Interpreter','latex')
hold off
    
lgd = legend([h1 h2 h3]);
lgd.NumColumns = 1;
lgd.Location = "northwest";



%%%%%%%% Time plot %%%%%%%%%%
fig = figure();
yyaxis right
h1 =  semilogx(tref,e11,'LineStyle', ':','Marker','none','LineWidth',2, ...
    'DisplayName',strcat('\epsilon_{input - experimental}'))  ;
hold all
trev = [tref+1e-10 ;flip(tref+1e-10)]';
inbetween = [(e11-e11Sd);flip(e11Sd +e11)];
h4 = fill(trev, inbetween, Orange,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
    'Marker','none',"DisplayName", "S. Deviation");
ylabel("$\epsilon$",'FontName','Times New Roman','FontSize', ...
    18,'Interpreter','latex')
xlim([1.9, max(tref)*1.2])
ylim([0 max(e11)*2])

yyaxis left
h2 = semilogx(tref,s11,'LineStyle', '-','Marker','none','LineWidth', ...
    1,'Color',Col(1),'DisplayName',strcat('\sigma_{experimental}'))  ;  
h3 =  semilogx(tref,s11_R,'LineStyle', ':','Marker','none', ...
    'LineWidth',2,'Color',Col(3),'DisplayName',strcat('\sigma_{model}'));
inbetween = [(s11-s11Sd);flip(s11Sd+s11)];
h4 = fill(trev, inbetween, Blue,'FaceAlpha',0.5, 'EdgeAlpha',0.2, ...
    'Marker','none',"DisplayName", "S. Deviation");
hold off
ylim([0 max(s11)*1.2])
set(gca,'fontname','Times New Roman')
fontsize(gca,16,"points")
set(fig,'defaultAxesColorOrder',[[0 0 1]; [1 1 1]]);
set(gcf, 'color', [1 1 1])
fig.Units = 'centimeters';
set(gcf,'Position',[3 3 19 15.2])
xlabel("t [s]",'FontName','Times New Roman','FontSize',18, ...
    'Interpreter','latex')
ylabel("$\sigma$ [MPa]",'FontName','Times New Roman','FontSize', ...
    18,'Interpreter','latex')
hold off
    
lgd = legend([h1 h2 h3]);
lgd.NumColumns = 1;
lgd.Location = "southeast";

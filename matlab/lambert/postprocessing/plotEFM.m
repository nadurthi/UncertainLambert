% plotEFM.m
% AUTHOR:           Robyn Woollands (robyn.woollands@gmail.com)
% DATE WRITTEN:     May, 2016
% LAST MODIFIED:    May, 2016
% AFFILIATION:      Department of Aerospace Engineering, Texas A & M University, College Station, TX
%
% DESCRIPTION:      Plots EFMs for number of solutions (possible/feasible)
%                   and delta V (total/departure,arrival)

close all
clear
clc

%% BEGIN User Input

possible_solns = 0;     % Plot EFM for number of possible solutions (YES = 1, NO = 0)
feasible_solns = 0;     % Plot EFM for number of feasible solutions (YES = 1, NO = 0)
delta_V        = 0;     % Plot EFM for total delta V (YES = 1, NO = 0)
delta_VD       = 0;     % Plot EFM for departure delta V (YES = 1, NO = 0)
delta_VA       = 0;     % Plot EFM for arrival delta V (YES = 1, NO = 0)
total_DV       = 1;     % Plot EFM for total minimum delta V (YES = 1, NO = 0)
total_DVD      = 1;     % Plot EFM for departure minimum delta V (YES = 1, NO = 0)
total_DVA      = 1;     % Plot EFM for arrival minimum delta V (YES = 1, NO = 0)

% END User Input

%% Load Data

load SOLN_p
DATAp = N_soln_p;

load SOLN_f
DATAf = N_soln_f;

load DelV
DATAv = DelV;

load DelVD
DATAvd = DelVD;

load DelVA
DATAva = DelVA;

load junkins_colormap.mat

%% Number of Possible Solutions
if possible_solns == 1
    SUMp = zeros(size(DATAp(:,:,1)));
    for i = 1:length(DATAp(1,1,:))
        SUMp = SUMp + DATAp(:,:,i);
    end
    saturate_level = max(max(max(SUMp)));
    
    for i = 1:length(DATAp(1,1,:))
        
        if i == 1
            XX = DATAp(:,:,i);
        elseif i > 1
            XX = XX + DATAp(:,:,i);
        end
        if mod(i,2) == 1
            figure('color',[1 1 1])
            set(gca, 'FontName', 'Helvetica','FontSize',16)
            imagesc(XX, [0, saturate_level] )
            set(gca,'YDir','normal')
            colorbar('FontName', 'Helvetica','FontSize',16)
            ylabel('Time of Flight (mins)')
            xlabel('t_0 (Mins Past Start)')
            title(['Number of Possible Solutions (N=',num2str(floor(i/2)),')'])
            saveas(gcf,['SolnPoss',num2str(ceil((i-1)/2)),'.png'])
            saveas(gcf,['SolnPoss',num2str(ceil((i-1)/2)),'.epsc'])
            saveas(gcf,['SolnPoss',num2str(ceil((i-1)/2)),'.fig'])
            saveas(gcf,['SolnPoss',num2str(ceil((i-1)/2)),'.jpg'])
        end
        
    end
end
%% Number of Feasible Solutions
if feasible_solns == 1
    SUMf = zeros(size(DATAf(:,:,1)));
    for i = 1:length(DATAf(1,1,:))
        SUMf = SUMf + DATAf(:,:,i);
    end
    saturate_level = max(max(max(SUMf)));
    
    for i = 1:length(DATAf(1,1,:))
        
        if i == 1
            XX = DATAf(:,:,i);
        elseif i > 1
            XX = XX + DATAf(:,:,i);
        end
        if mod(i,2) == 1
            figure('color',[1 1 1])
            imagesc(XX, [0, saturate_level] )
            set(gca, 'FontName', 'Helvetica','FontSize',16)
            set(gca,'YDir','normal')
            colorbar('FontName', 'Helvetica','FontSize',16)
            ylabel('Time of Flight (mins)')
            xlabel('t_0 (Mins Past Start)')
            title(['Number of Feasible Solutions (N=',num2str(i/2),')'])
            
            saveas(gcf,['SolnFeas',num2str(ceil((i-1)/2)),'.png'])
            saveas(gcf,['SolnFeas',num2str(ceil((i-1)/2)),'.epsc'])
            saveas(gcf,['SolnFeas',num2str(ceil((i-1)/2)),'.fig'])
            saveas(gcf,['SolnFeas',num2str(ceil((i-1)/2)),'.jpg'])
        end
        
    end
end
%% Total Delta V
saturate_level = max(max(max(DATAv)));
if delta_V == 1
    
    for i = 1:length(DATAv(1,1,:))
        
        figure('color',[1 1 1])
        imagesc(DATAv(:,:,i), [0, saturate_level] )
        set( gcf, 'Colormap', junkins_colormap);
        set(gca,'YDir','normal','FontName', 'Helvetica','FontSize',16)
        colorbar
        ylabel('Time of Flight (mins)')
        xlabel('t_0 (Mins Past Start)')
        
        if i == 1
            title(['Total \Delta V (N=',num2str(ceil((i-1)/2)),')'])
            saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'.png'])
            saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'.epsc'])
            saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'.fig'])
            saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'.jpg'])
        end
        if i > 1
            if mod(i,2) == 0
                title(['Total \Delta V (N=',num2str(ceil((i-1)/2)),', Upper Branch)'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.png'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.epsc'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.fig'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.jpg'])
            end
            if mod(i,2) == 1
                title(['Total \Delta V (N=',num2str(ceil((i-1)/2)),', Lower Branch)'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.png'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.fig'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.jpg'])
            end
        end
    end
    
end

if total_DV == 1
    
    % Total Delta V
    EFM = zeros(size(DATAv,1),size(DATAv,2));
    MIN = saturate_level;
    for i = 1:size(DATAv,1)
        for j = 1:size(DATAv,2)
            for k = 1:size(DATAv,3)
                if DATAv(i,j,k) < MIN && DATAv(i,j,k) > 0
                    EFM(i,j) = DATAv(i,j,k);
                    MIN      = EFM(i,j);
                end
            end
            MIN = saturate_level;
        end
    end
    
    figure('color',[1 1 1])
    imagesc(EFM, [0, saturate_level] )
    set( gcf, 'Colormap', junkins_colormap);
    set(gca,'YDir','normal','FontName', 'Helvetica','FontSize',16)
    colorbar
    ylabel('Time of Flight (mins)')
    xlabel('t_0 (Mins Past Start)')
    
    title('Total \Delta V')
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.png'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.fig'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.jpg'])

end

%% Departure Delta V
saturate_level = max(max(max(DATAvd)));
if delta_VD == 1
    
    for i = 1:length(DATAvd(1,1,:))
        
        figure('color',[1 1 1])
        set(gca, 'FontName', 'Helvetica','FontSize',16)
        imagesc(DATAvd(:,:,i), [0, saturate_level] )
        set( gcf, 'Colormap', junkins_colormap);
        set(gca,'YDir','normal')
        colorbar('FontName', 'Helvetica','FontSize',16)
        ylabel('Time of Flight (mins)')
        xlabel('t_0 (Mins Past Start)')
        if i == 1
            title(['Departure \DeltaV (N=',num2str(ceil((i-1)/2)),')'])
            saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'.png'])
            saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'.epsc'])
            saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'.fig'])
            saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'.jpg'])
        end
        if i > 1
            if mod(i,2) == 0
                title(['Departure \DeltaV (N=',num2str(ceil((i-1)/2)),', Upper Branch)'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_upper.png'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_upper.epsc'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.fig'])
                saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_upper.jpg'])
            end
            if mod(i,2) == 1
                title(['Departure \Delta V (N=',num2str(ceil((i-1)/2)),', Lower Branch)'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_lower.png'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_lower.fig'])
                saveas(gcf,['DVD_N',num2str(ceil((i-1)/2)),'_lower.jpg'])
            end
        end
        
    end

    % Total Delta V
    EFM = zeros(size(DATAvd,1),size(DATAvd,2));
    MIN = saturate_level;
    for i = 1:size(DATAvd,1)
        for j = 1:size(DATAvd,2)
            for k = 1:size(DATAvd,3)
                if DATAvd(i,j,k) < MIN && DATAvd(i,j,k) > 0
                    EFM(i,j) = DATAvd(i,j,k);
                    MIN      = EFM(i,j);
                end
            end
            MIN = saturate_level;
        end
    end
    
end

if total_DVD == 1
    
    figure('color',[1 1 1])
    imagesc(EFM, [0, saturate_level] )
    set( gcf, 'Colormap', junkins_colormap);
    set(gca,'YDir','normal','FontName', 'Helvetica','FontSize',16)
    colorbar
    ylabel('Time of Flight (mins)')
    xlabel('t_0 (Mins Past Start)')
    
    title('Departure \Delta V')
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.png'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.fig'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.jpg'])
    
end

%% Arrival Delta
saturate_level = max(max(max(DATAva)));
if delta_VA == 1
    
    for i = 1:length(DATAva(1,1,:))
        
        figure('color',[1 1 1])
        imagesc(DATAva(:,:,i), [0, saturate_level] )
        set( gcf, 'Colormap', junkins_colormap);
        set(gca,'YDir','normal')
        colorbar
        ylabel('Time of Flight (mins)')
        xlabel('t_0 (Mins Past Start)')
        if i == 1
            title(['Arrival \Delta V (N=',num2str(ceil((i-1)/2)),')'])
            saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'.png'])
            saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'.epsc'])
            saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'.fig'])
            saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'.jpg'])
        end
        if i > 1
            if mod(i,2) == 0
                title(['Arrival \Delta V (N=',num2str(ceil((i-1)/2)),', Upper Branch)'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_upper.png'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_upper.epsc'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_upper.fig'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_upper.jpg'])
            end
            if mod(i,2) == 1
                title(['Arrival \Delta V (N=',num2str(ceil((i-1)/2)),', Lower Branch)'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_lower.png'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_lower.fig'])
                saveas(gcf,['DVA_N',num2str(ceil((i-1)/2)),'_lower.jpg'])
            end
        end
        
    end
    
end

if total_DVA == 1
    
    % Total Delta V
    EFM = zeros(size(DATAva,1),size(DATAva,2));
    MIN = saturate_level;
    for i = 1:size(DATAva,1)
        for j = 1:size(DATAva,2)
            for k = 1:size(DATAva,3)
                if DATAva(i,j,k) < MIN && DATAva(i,j,k) > 0
                    EFM(i,j) = DATAva(i,j,k);
                    MIN      = EFM(i,j);
                end
            end
            MIN = saturate_level;
        end
    end
    
    figure('color',[1 1 1])
    imagesc(EFM, [0, saturate_level] )
    set( gcf, 'Colormap', junkins_colormap);
    set(gca,'YDir','normal','FontName', 'Helvetica','FontSize',16)
    colorbar
    ylabel('Time of Flight (mins)')
    xlabel('t_0 (Mins Past Start)')
    
    title('Arrival \Delta V')
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.png'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.epsc'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.fig'])
    saveas(gcf,['DVT_N',num2str(ceil((i-1)/2)),'_lower.jpg'])
    
end
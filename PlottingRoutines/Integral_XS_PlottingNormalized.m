cd('~/Documents/Juno/Jup_Hyd/XS/');
% fname = 'OG_Integral_XS_CTMC.txt';
fname = 'OG_Integral_XS_Normalized.txt';
% fname1 = 'OG_Integral_XS_CTMC_Interpolated.dat';
fname1 = 'Integral_XS_Normalized_Interpolated.dat';
fid = fopen(fname);
fid1 = fopen(fname1);

nProc = 0;
nProcN = 0;

% if contains(fname,'CTMC')
    while ~feof(fid)
        nEng = 0;
        tline = fgetl(fid);
        if contains(tline,'H^{')
            nProc = nProc + 1;
            proc = regexp(tline,'\s','split');
            procNames(nProc) = string(proc{1});
            fgetl(fid);
            fgetl(fid);
            tline = fgetl(fid);
        end
        while ~contains(tline,'!') && ~feof(fid)
            nEng = nEng + 1;
            tmp = strsplit(tline);
            while isempty(tmp{1})
                tmp = tmp(2:end);
            end
            energy(nEng,nProc) = str2double(tmp{1});
            xs(nEng,1,nProc) = str2double(tmp{2});
            xs(nEng,2,nProc) = str2double(tmp{4});
            xs(nEng,3,nProc) = str2double(tmp{6});
            tline = fgetl(fid);
        end
    end
% else
    while ~feof(fid1)
        nEng = 0;
        tline = fgetl(fid1);
        if contains(tline,'H^{')
            nProcN = nProcN + 1;
            proc = regexp(tline,'\s','split');
            procNames(nProcN) = string(proc{1});
            fgetl(fid1);
            fgetl(fid1);
            tline = fgetl(fid1);
        end
        while ~contains(tline,'!') && ~feof(fid1)
            nEng = nEng + 1;
            tmp = strsplit(tline);
            while isempty(tmp{1})
                tmp = tmp(2:end);
            end
            energyN(nEng,nProcN) = str2double(tmp{1});
            xsN(nEng,1,nProcN) = str2double(tmp{2});
            xsN(nEng,2,nProcN) = str2double(tmp{3});
            xsN(nEng,3,nProcN) = str2double(tmp{4});
            tline = fgetl(fid1);
        end
    end
% end
fclose(fid);
fclose(fid1);

% figure(1)
% colors = [1 0 0; 0 1 0; 0 0 1]; 
% for i = 1:nProc
%     subplot(2,5,i)
%     hold on
%     for j = 1:3
%         loglog(energy(:,j),xs(:,j,i),'Color',colors(j,:),'LineWidth',2,'MarkerSize',6)
%     end
%     plot(energyN(:,1),xsN(:,1,i),'k','LineWidth',2,'LineStyle','--')
%     for j = 1:3
%         loglog(energyN(:,j),xsN(:,j,i),'Color',colors(j,:),'LineWidth',2,'LineStyle','--')
%     end
%     title(procNames{i});
% %     if i == 1
%         legend({'H^{+}','H','H^{-}','Interpolated'})
% %     end
%     xlabel('Projectile Energy [keV]')
%     ylabel('Cross-Section [cm^{-2}]')
%     xticks([1 10 100 1000 10000 100000])
% %     ylim([max(max(xs(:,:,i)))/1000 max(max(xs(:,:,i)))])
%     xlim([1 1e5])
%     set(gca,'FontSize',16,'FontWeight','bold','XMinorTick','on','XScale',...
%     'log','YMinorTick','on','YScale','log')
% %     axis square
%     hold off
% %     keyboard
% end
% 
% figure(2)
% colors = [0 0 0;colormap(jet(9))];
% titles = ["H^{-}" "H" "H^{+}"];
% for i = 1:3
%     subplot(1,3,i)
%     hold on
%     for j = 1:nProc
%         loglog(energyN(:,i),xsN(:,i,j),'Color',colors(j,:),'LineWidth',2)
%     end
%     legend(procNames)
%     title(titles{i})
%     xlabel('Projectile Energy [keV]')
%     ylabel('Cross-Section [cm^{-2}]')
%     xticks([1 10 100 1000 10000 100000])
%     ylim([1e-25 1e-14])
%     set(gca,'FontSize',16,'FontWeight','bold','XMinorTick','on','XScale',...
%     'log','YMinorTick','on','YScale','log')
%     hold off
% end

ct = xsN(:,:,3)+xsN(:,:,6)+xsN(:,:,7);
strip = xsN(:,:,4)+xsN(:,:,5);
other = sum(xsN,3) - ct - strip;

figure(3)
colors = [0 0 0;1 0 0;0 0 1];
titles = ["H^{-}" "H" "H^{+}"];
for i = 1:3
    subplot(1,3,i)
    hold on
    
    loglog(energyN(:,i),other(:,i),'Color',colors(1,:),'LineWidth',2)
    loglog(energyN(:,i),ct(:,i),'Color',colors(2,:),'LineWidth',2)
    loglog(energyN(:,i),strip(:,i),'Color',colors(3,:),'LineWidth',2)

    legend({'Other - SI + DI + TEX + PEX + ES','Charge Transfer - TI + SC + DC','Stripping - SS + DS'})
    title(titles{i})
    xlabel('Projectile Energy [keV]')
    ylabel('Cross-Section [cm^{-2}]')
    xticks([1 10 100 1000 10000 100000])
    ylim([1e-25 1e-14])
    set(gca,'FontSize',16,'FontWeight','bold','XMinorTick','on','XScale',...
    'log','YMinorTick','on','YScale','log')
    hold off
end

% fgetl(fid)
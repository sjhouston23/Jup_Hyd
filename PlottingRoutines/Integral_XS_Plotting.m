cd('~/Documents/Juno/Jup_Hyd/XS/');
fname = 'OG_Integral_XS_CTMC.txt';
fname1 = 'OG_Integral_XS_Normalized.txt';
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
            xs(nEng,2,nProc) = str2double(tmp{3});
            xs(nEng,3,nProc) = str2double(tmp{4});
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
            xsN(nEng,2,nProcN) = str2double(tmp{4});
            xsN(nEng,3,nProcN) = str2double(tmp{6});
            tline = fgetl(fid1);
        end
    end
% end
fclose(fid);
fclose(fid1);

figure(1)
colors = [1 0 0; 0 1 0; 0 0 1]; 
for i = 1:nProc
    subplot(2,5,i)
    hold on
    for j = 1:3
        loglog(energy(:,j),xs(:,j,i),'Color',colors(j,:),'LineWidth',2,'MarkerSize',6)
    end
    plot(energyN(:,1),xsN(:,1,i),'k','LineWidth',2,'LineStyle','--')
    for j = 1:3
        loglog(energyN(:,j),xsN(:,j,i),'Color',colors(j,:),'LineWidth',2,'LineStyle','--')
    end
    title(procNames{i});
%     if i == 1
        legend({'H^{+}','H','H^{-}','Normalized'})
%     end
    xlabel('Projectile Energy [keV]')
    ylabel('Cross-Section [cm^{-2}]')
    xticks([1 10 100 1000 10000 100000])
%     ylim([1e-20 1e-14])
    xlim([1 1e5])
    set(gca,'FontSize',16,'FontWeight','bold','XMinorTick','on','XScale',...
    'log','YMinorTick','on','YScale','log'),
%     axis square
    hold off
end

% fgetl(fid)
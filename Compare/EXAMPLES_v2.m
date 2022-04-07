%% Two-dimensional X-data
%clear, clc, close all
%{
load spectra
tt = repmat((1:3)',20,1);
ncomp = 10;
Xcal = NIR(tt<3,:);  Ycal = octane(tt<3,:);
Xval = NIR(tt==3,:); Yval = octane(tt==3,:);
Ysec = log(abs(Ycal));
[B, T, W, Q, R, Wa, Pt] = ncpls(ncomp, Xcal, Ycal, Ysec);
[B1, W1, P1, T1] = cpls(ncomp, Xcal, Ycal, Ysec);
[B2, W2, P2, T2] = cpls0(ncomp, Xcal, Ycal, Ysec);
ypred = (Xval - mean(Xcal,1))*B + mean(Ycal);
RMSEP = zeros(ncomp,size(Yval,2)); R2 = RMSEP;
for i=1:ncomp
    RMSEP(i,:) = sqrt(mean((Yval-ypred(:,i)).^2));
    R2(i,:) = (1 - mean((Yval-ypred(:,i)).^2)./var(Yval))*100;
end
figure
plot(R2(:,1), '-r')
ylabel('R^2')
xlabel('#comp')
grid on
%}

% {
clear, clc, close all
load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Multiway\sugar.mat')
Yprim = y(:,1:2); Ysec = y(:,3:end);
% load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Dough.mat')
% Yprim = Ytrain(:,1); Ysec = Ytrain(:,2:end); X = Xtrain;

pc = 5;
tic, [B, W, P, T, mX, cc]       = cpls(pc, X, Yprim, Ysec); tid1 = toc
tic, [B0, W0, P0, T0, mX0, cc0] = cpls0(pc, X, Yprim, Ysec); tid0 = toc

%figure, plot(B(2:end,:)), grid on
%figure, plot(B0(2:end,:)), grid on
%norm(B(2:end,:)-B0(2:end,:),'fro')
figure, plot(B(2:end,:,1)), grid on
figure, plot(B0(2:end,:,1)), grid on
norm(B(2:end,:,1)-B0(2:end,:,1),'fro')

%}


%% Effect of orthogonalization of w_j and unfolded analysis
% Three-dimensional X-data
% clear, clc, close all
if ~isempty(dir('C:\Users\kristl'))
    load('C:\Users\kristl\Dropbox\Datasett\Multiway\sugar.mat')
elseif ~isempty(dir('C:\Users\ulfin'))
    load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Multiway\sugar.mat')
else
    load('~/Dropbox/Datasett/Multiway/sugar.mat')
end
X = reshape(X,DimX);

if 1==2 % Remove Rayleigh scattering
    X(:,100:180,1) = 0;
    X(:,70:150,2) = 0;
    X(:,30:120,3) = 0;
    X(:,1:90,4) = 0;
    X(:,1:60,5) = 0;
end

ncomp = 25;

Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,2:3); YcalAdd = [];%y(1:2:end,1);
Xval = X(2:2:end,:,:); Yval = y(2:2:end,2:3);

% Orth w_j, folded analysis
[Ba,Ta,Wa,Qa,Ra,Woa,Pta] = ncpls(ncomp, Xcal, Ycal, YcalAdd);
ypreda = GMP(Xval-mean(Xcal,1), Ba, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% Non-orth w_j, folded analysis
[Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false);
ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% Orth w_j, unfolded analysis
[Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% % Non-orth w_j, unfolded analysis
% [Bd,Td,Wd,Qd,Rd,Wod] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false, false);
% ypredd = GMP(Xval-mean(Xcal,1), Bd, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);

RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
for i=1:ncomp
    RMSEPa(i,:) = sqrt(mean((Yval-squeeze(ypreda(:,i,:))).^2));
    R2a(i,:) = (1 - mean((Yval-squeeze(ypreda(:,i,:))).^2)./var(Yval))*100;
    RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
    RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
    R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
end

figure
c = colormap('lines');
hold on
plot(R2c(:,1), '-', 'Color',c(1,:),'LineWidth',2)
plot(R2b(:,1), '--', 'Color',c(1,:),'LineWidth',2)
plot(R2a(:,1), ':','Color',c(1,:),'LineWidth',2)
plot(R2c(:,2), '-','Color',c(2,:),'LineWidth',2)
plot(R2b(:,2), '--','Color',c(2,:),'LineWidth',2)
plot(R2a(:,2), ':','Color',c(2,:),'LineWidth',2)
grid on
% title('Sugar')
legend('Colour - unfolded', 'Colour - multilinear', 'Colour - multilinear+orthogonal', 'Ash - unfolded', 'Ash - multilinear', 'Ash - multilinear+orthogonal', 'Location', 'SouthEast')
xlabel('# comp.')
ylabel('% explained variance')
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])

figure
P = get(gcf,'Position'); set(gcf,'Position', [P(1:2) 560, 145])
subplot(133)
imagesc(Wa{1}'*Wa{1})
title("Mul.+orth. w^{[1]} : w^{[1]}'w^{[1]}")
subplot(132)
imagesc(Wb{1}'*Wb{1})
title("Multilinear w^{[1]} : w^{[1]}'w^{[1]}")
subplot(131)
vec = @(x)x(:);
ww = [];
for i=1:ncomp
    ww = [ww vec(Woc(:,:,i))];
end
imagesc(ww'*ww)
title("Unfolded W : W'W")
colorbar
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])

for i=1:ncomp
    Wodd1 = GOuter({Wc{1}(:,i),Wc{2}(:,i)});
    Wodd2 = Woc(:,:,i);
    disp(corr(Wodd1(:),Wodd2(:)))
    disp(norm(Wodd2-Wodd1, 'fro'))
end

fig = figure;
fig.Position = [264   446   780   420];
subplot(2,3,4:5)
plot(EmAx,Wb{1}(:,1:2), 'LineWidth',2)
axis tight
xlabel('emission (nm)')
ylabel('intensity')
title('w^{[2]}')
legend('Component 1', 'Component 2')
subplot(2,3,1:2)
plot(ExAx,Wb{2}(:,1:2), 'LineWidth',2)
axis tight
xlabel('excitation (nm)')
ylabel('intensity')
title('w^{[1]}')
legend('Component 1', 'Component 2', 'Location','SouthWest')
subplot(2,3,3)
mesh(-Wob(:,:,1))
ylabel('emission (nm)')
xlabel('exitation (nm)')
zlabel('intensity')
title('W^O, component 1')
subplot(2,3,6)
mesh(-Wob(:,:,2))
ylabel('emission (nm)')
xlabel('excitation (nm)')
zlabel('intensity')
title('W^O, component 2')
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])


%% NPLS, NCPLS, NCPLS-add  -  sugar
if ~isempty(dir('C:\Users\kristl'))
    load('C:\Users\kristl\Dropbox\Datasett\Multiway\sugar.mat')
    addpath(genpath('..\NPLS')) % Forutsetter at arbeidskatalogen er samme som dette skriptet
elseif ~isempty(dir('C:\Users\ulfin'))
    load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Multiway\sugar.mat')
    addpath(genpath('..\NPLS')) % Forutsetter at arbeidskatalogen er samme som dette skriptet
else
    load('~/Dropbox/Datasett/Multiway/sugar.mat')
    addpath(genpath('../NPLS')) % Forutsetter at arbeidskatalogen er samme som dette skriptet
end
X = reshape(X,DimX);

ncomp = 20;

t_str = num2str(time);
month = dummyvar(str2num(t_str(:,1)));
year = dummyvar(str2num(t_str(:,2:3)));
period = dummyvar(str2num(t_str(:,4)));

% 'time gives the time of sampling xyyz. x is the month (2:Oct,3:Nov,4:Dec,5:Jan),yy the day, and z the time a day (1:Morning,2:Afternoon,3:Nigth)'
Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,3); YcalAdd = [y(1:2:end,2), year(1:2:end,:)];
Xval = X(2:2:end,:,:); Yval = y(2:2:end,3);
%Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,2); YcalAdd = y(1:2:end,[1,3]);
%Xval = X(2:2:end,:,:); Yval = y(2:2:end,2);
%Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,1); YcalAdd = y(1:2:end,2:3);
%Xval = X(2:2:end,:,:); Yval = y(2:2:end,1);

n_time = 100;
time_sugar = zeros(n_time,4);
for t=1:n_time % timing
    % NPLS
    tic
    [Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
    RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
    for i=1:ncomp
        [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
        RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
        R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
    end
    time_sugar(t,1) = toc;
    % NCPLS
    tic
    [Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, [], true, false);
    ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
    for i=1:ncomp
        RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
        R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
    end
    time_sugar(t,2) = toc;
    % NCPLS - additional
    tic
    [Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
    ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
    for i=1:ncomp
        RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
        R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
    end
    time_sugar(t,3) = toc;
    % NCPLS - multilinear
    tic
    [Bd,Td,Wd,Qd,Rd,Wod] = ncpls(ncomp, Xcal, Ycal, [], true, true);
    ypredd = GMP(Xval-mean(Xcal,1), Bd, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPd = zeros(ncomp,size(Yval,2)); R2d = RMSEPd;
    for i=1:ncomp
        RMSEPd(i,:) = sqrt(mean((Yval-squeeze(ypredd(:,i,:))).^2));
        R2d(i,:) = (1 - mean((Yval-squeeze(ypredd(:,i,:))).^2)./var(Yval))*100;
    end
    time_sugar(t,4) = toc;
end

% Explained variance
figure, hold on
c = colormap('lines');
plot(R2a(:,1), '-', 'Color',c(1,:),'LineWidth',2);% plot(R2a(:,2), ']--k')
plot(R2d(:,1), '--', 'Color',c(1,:),'LineWidth',2);% plot(R2b(:,2), '--r')
plot(R2b(:,1), ':', 'Color',c(1,:),'LineWidth',2);% plot(R2c(:,2), '--b')
plot(R2c(:,1), ':', 'Color',c(2,:),'LineWidth',2);% plot(R2c(:,2), '--b')
xlabel('# comp.')
ylabel('% explained variance')
grid on
legend('N-PLS','N-CPLS - multilinear','N-CPLS - unfolded','N-CPLS - unfolded with Y_{additional}','Location','SouthEast')
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])
% legend('NPLS, Resp. 1', 'NPLS, Resp. 2', 'NCPLS, Resp. 1', 'NCPLS, Resp. 2', 'NCPLS-add, Resp. 1', 'NCPLS-add, Resp. 2','Location','SouthEast')


%% NPLS, NCPLS, NCPLS-add  -  milk fats, single response
if ~isempty(dir('C:\Users'))
    load('..\Mishra\Problem 3 Multiway CPLS\ncpls_PM\milk_fats.mat')
else
    load('../Mishra/Problem 3 Multiway CPLS/ncpls_PM/milk_fats.mat')
end
X = X{2};

ncomp = 10;

Xcal = X(1:2:end,:,:); Ycal = fats(1:2:end,1); YcalAdd = [protein(1:2:end,1),total_solids(1:2:end,1)];
Xval = X(2:2:end,:,:); Yval = fats(2:2:end,1);

n_time = 100;
time_milk = zeros(n_time,4);
for t=1:n_time % timing
    % NPLS
    tic
    [Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
    RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
    for i=1:ncomp
        [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
        RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
        R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
    end
    time_milk(t,1) = toc;
    % NCPLS
    tic
    [Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, [], true, false);
    ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
    for i=1:ncomp
        RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
        R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
    end
    time_milk(t,2) = toc;
    % NCPLS - additional
    tic
    [Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
    ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
    for i=1:ncomp
        RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
        R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
    end
    time_milk(t,3) = toc;
    % NCPLS - multilinear
    tic
    [Bd,Td,Wd,Qd,Rd,Wod] = ncpls(ncomp, Xcal, Ycal, [], true, true);
    ypredd = GMP(Xval-mean(Xcal,1), Bd, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPd = zeros(ncomp,size(Yval,2)); R2d = RMSEPd;
    for i=1:ncomp
        RMSEPd(i,:) = sqrt(mean((Yval-squeeze(ypredd(:,i,:))).^2));
        R2d(i,:) = (1 - mean((Yval-squeeze(ypredd(:,i,:))).^2)./var(Yval))*100;
    end
    time_milk(t,4) = toc;
end

% Explained variance
figure, hold on
c = colormap('lines');
plot(R2a(:,1), '-', 'Color',c(1,:),'LineWidth',2);% plot(R2a(:,2), ']--k')
plot(R2d(:,1), '--', 'Color',c(1,:),'LineWidth',2);% plot(R2b(:,2), '--r')
plot(R2b(:,1), ':', 'Color',c(1,:),'LineWidth',2);% plot(R2c(:,2), '--b')
plot(R2c(:,1), ':', 'Color',c(2,:),'LineWidth',2);% plot(R2c(:,2), '--b')
xlabel('# comp.')
ylabel('% explained variance')
grid on
legend('N-PLS','N-CPLS - multilinear','N-CPLS - unfolded','N-CPLS - unfolded with Y_{additional}','Location','SouthEast')
ylim([0,100])
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])
% legend('NPLS, Resp. 1', 'NPLS, Resp. 2', 'NCPLS, Resp. 1', 'NCPLS, Resp. 2', 'NCPLS-add, Resp. 1', 'NCPLS-add, Resp. 2','Location','SouthEast')



% %% NPLS, NCPLS, NCPLS-add  -  milk fats, double response
% if ~isempty(dir('C:\Users'))
%     load('..\Mishra\Problem 3 Multiway CPLS\ncpls_PM\milk_fats.mat')
% else
%     load('../Mishra/Problem 3 Multiway CPLS/ncpls_PM/milk_fats.mat')
% end
% X = X{2};
% 
% ncomp = 25;
% 
% Xcal = X(1:2:end,:,:); Ycal = [fats(1:2:end,1), protein(1:2:end,1)]; YcalAdd = total_solids(1:2:end,1);
% Xval = X(2:2:end,:,:); Yval = [fats(2:2:end,1), protein(1:2:end,1)];
% 
% % NPLS
% [Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
% % NCPLS
% [Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, [], true, false);
% ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% % NCPLS - additional
% [Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
% ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% 
% RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
% RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
% RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
% for i=1:ncomp
%     [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
%     RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
%     R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
%     RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
%     R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
%     RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
%     R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
% end
% 
% % Explained variance
% figure, hold on
% plot(R2a(:,1), '-k'),% plot(R2a(:,2), '--k')
% plot(R2b(:,1), '-r'),% plot(R2b(:,2), '--r')
% plot(R2c(:,1), '-b'),% plot(R2c(:,2), '--b')
% xlabel('# comp.')
% ylabel('% explained variance')
% grid on
% legend('Multiway PLS','Multiway CPLS','Multiway CPLS (Yadd)','Location','SouthEast')
% ylim([0,100])
% pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])
% % legend('NPLS, Resp. 1', 'NPLS, Resp. 2', 'NCPLS, Resp. 1', 'NCPLS, Resp. 2', 'NCPLS-add, Resp. 1', 'NCPLS-add, Resp. 2','Location','SouthEast')



% %% NPLS, NCPLS, NCPLS-add  -  dorrit
% load('dorrit.mat')
% X = EEM.data;
% Y = Y.data;
% 
% ncomp = 25;
% 
% Xcal = X(1:2:end,:,:); Ycal = Y(1:2:end,1); YcalAdd = Y(1:2:end,2:end);
% Xval = X(2:2:end,:,:); Yval = Y(2:2:end,1);
% 
% % NPLS
% [Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
% % NCPLS
% [Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, [], true, false);
% ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% % NCPLS - additional
% [Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
% ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% 
% RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
% RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
% RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
% for i=1:ncomp
%     [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
%     RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
%     R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
%     RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
%     R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
%     RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
%     R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
% end
% 
% % Explained variance
% figure, hold on
% plot(R2a(:,1), '-k'),% plot(R2a(:,2), '--k')
% plot(R2b(:,1), '-r'),% plot(R2b(:,2), '--r')
% plot(R2c(:,1), '-b'),% plot(R2c(:,2), '--b')
% xlabel('# comp.')
% ylabel('% explained variance')
% grid on
% legend('Multiway PLS','Multiway CPLS','Multiway CPLS (Yadd)','Location','SouthEast')
% ylim([0,100])
% pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])
% % legend('NPLS, Resp. 1', 'NPLS, Resp. 2', 'NCPLS, Resp. 1', 'NCPLS, Resp. 2', 'NCPLS-add, Resp. 1', 'NCPLS-add, Resp. 2','Location','SouthEast')




%% NPLS, NCPLS, NCPLS-add  -  ORL faces
load('orl.mat')
X = reshape(Xtrain, [400,112,92]);
Y = dummyvar(Ytrain);

ncomp = 50;

Xcal = X(1:2:end,:,:); Ycal = Y(1:2:end,:);
Xval = X(2:2:end,:,:); Yval = Y(2:2:end,:);
[u,s,v] = svds(Xtrain(1:2:end,:)-mean(Xtrain(1:2:end,:)),10);
YcalAdd = u*s;


% Effect of restrictions
figure
P = get(gcf,'Position'); set(gcf,'Position', [P(1:2) 400, 400])
subplot(221)
imagesc(squeeze(Xcal(1,:,:)));colormap(bone)
subplot(222)
imagesc(squeeze(Xcal(10,:,:)));colormap(bone)
subplot(223)
[Bb,Tb,Wb,Qb,Rb,Wob,Ptb] = ncpls(2, Xcal, Ycal, [], true, false);
imagesc(squeeze(Wob(:,:,1)));colormap(bone)
title('Unfolded')
subplot(224)
[Bb,Tb,Wb,Qb,Rb,Wob,Ptb] = ncpls(2, Xcal, Ycal, [], false, true);
imagesc(squeeze(Wob(:,:,1)));colormap(bone)
title('Multilinear')
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])

n_time = 10;
time_ORL = zeros(n_time,4);
for t=1:n_time % timing
    % NPLS
    tic
    [Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
    RMSEPa = zeros(ncomp,1); R2a = RMSEPa;
    for i=1:ncomp
        [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
        [~,ma] = max(squeeze(ypreda),[],2);
        [~,m2] = max(Yval,[],2);
        R2a(i,:) = sum(ma==m2)/size(Yval,1)*100;%sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    end
    time_ORL(t,1) = toc;
    % NCPLS
    tic
    [Bb,Tb,Wb,Qb,Rb,Wob,Ptb] = ncpls(ncomp, Xcal, Ycal, [], true, false);
    ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPb = zeros(ncomp,1); R2b = RMSEPb;
    for i=1:ncomp
        [~,mb] = max(squeeze(ypredb(:,i,:)),[],2);
        [~,m2] = max(Yval,[],2);
        R2b(i,:) = sum(mb==m2)/size(Yval,1)*100;%sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    end
    time_ORL(t,2) = toc;
    % NCPLS - additional
    tic
    [Bc,Tc,Wc,Qc,Rc,Woc,Ptc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
    ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPc = zeros(ncomp,1); R2c = RMSEPc;
    for i=1:ncomp
        [~,mc] = max(squeeze(ypredc(:,i,:)),[],2);
        [~,m2] = max(Yval,[],2);
        R2c(i,:) = sum(mc==m2)/size(Yval,1)*100;%sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    end
    time_ORL(t,3) = toc;
    % NCPLS
    tic
    [Bd,Td,Wd,Qd,Rd,Wod,Ptd] = ncpls(ncomp, Xcal, Ycal, [], false, true);
    ypredd = GMP(Xval-mean(Xcal,1), Bd, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
    RMSEPd = zeros(ncomp,1); R2d = RMSEPd;
    for i=1:ncomp
        [~,md] = max(squeeze(ypredd(:,i,:)),[],2);
        [~,m2] = max(Yval,[],2);
        R2d(i,:) = sum(md==m2)/size(Yval,1)*100;%sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    end
    time_ORL(t,4) = toc;
end

% Explained variance
figure, hold on
c = colormap('lines');
plot(R2a(:,1), '-', 'Color',c(1,:),'LineWidth',2);% plot(R2a(:,2), ']--k')
plot(R2d(:,1), '--', 'Color',c(1,:),'LineWidth',2);% plot(R2b(:,2), '--r')
plot(R2b(:,1), ':', 'Color',c(1,:),'LineWidth',2);% plot(R2c(:,2), '--b')
plot(R2c(:,1), ':', 'Color',c(2,:),'LineWidth',2);% plot(R2c(:,2), '--b')
xlabel('# comp.')
ylabel('% correctly classified')
grid on
legend('N-PLS','N-CPLS - multilinear','N-CPLS - unfolded','N-CPLS - unfolded with Y_{additional}','Location','SouthEast','AutoUpdate',false)
i81 = find(R2a(:,1)>=80,1); i82 = find(R2b(:,1)>=80,1); i83 = find(R2c(:,1)>=80,1); i84 = find(R2d(:,1)>=80,1);
i91 = find(R2a(:,1)>=90,1); i92 = find(R2b(:,1)>=90,1); i93 = find(R2c(:,1)>=90,1); i94 = find(R2d(:,1)>=90,1);
plot([i81,i82,i83,i84,i91,i92,i93,i94], [R2a(i81,1),R2b(i82,1),R2c(i83,1),R2d(i84,1),R2a(i91,1),R2b(i92,1),R2c(i93,1),R2d(i94,1)],'xk')
plot([i81,i82,i83,i84,i91,i92,i93,i94], [R2a(i81,1),R2b(i82,1),R2c(i83,1),R2d(i84,1),R2a(i91,1),R2b(i92,1),R2c(i93,1),R2d(i94,1)],'or')
% legend('Multiway PLS','Multiway CPLS','Location','SouthEast')
ylim([0,100])
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])




%% Raman milk, various preps.
if(~isempty(dir('C:\Users\kristl')))
    load('C:\Users\kristl\Dropbox\Datasett\Ramanmilk_cleaned.mat')
else
    load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Ramanmilk_cleaned.mat')
end

Raman = Raman(:,shift>=120 & shift <=3100);
shift = shift(shift>=120 & shift <=3100);
X = zeros([size(Raman),9]);
X(:,:,1) = Raman;
for t = 1:8
    [X(:,:,t+1), parameters, par_names, model] = emsc(Raman, 'reference',Raman(579,:), 'terms',t);
end

ncomp = 25;

Xcal = X(cvseg<163,:,:); Ycal = CLA(cvseg<163,1); YcalAdd = Iodine(cvseg<163,1);
Xval = X(cvseg>=163,:,:); Yval = CLA(cvseg>=163,1);

% NPLS
[Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
% NCPLS
[Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, [], true, false);
ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% NCPLS - additional
[Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% NCPLS - additional - single method
[Bd,Td,Wd,Qd,Rd,Wod] = ncpls(ncomp, Xcal(:,:,7), Ycal, YcalAdd, true, false);
ypredd = GMP(Xval(:,:,7)-mean(Xcal(:,:,7),1), Bd, 1)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);

RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
RMSEPd = zeros(ncomp,size(Yval,2)); R2d = RMSEPd;
for i=1:ncomp
    [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
    RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
    R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
    RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
    RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
    R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
    RMSEPd(i,:) = sqrt(mean((Yval-squeeze(ypredd(:,i,:))).^2));
    R2d(i,:) = (1 - mean((Yval-squeeze(ypredd(:,i,:))).^2)./var(Yval))*100;
end

% Explained variance
figure, hold on
plot(R2a(:,1), '-k'),% plot(R2a(:,2), '--k')
plot(R2b(:,1), '-r'),% plot(R2b(:,2), '--r')
plot(R2c(:,1), '-b'),% plot(R2c(:,2), '--b')
plot(R2d(:,1), '--k'),% plot(R2c(:,2), '--b')
xlabel('# comp.')
ylabel('% explained variance')
grid on
legend('Multiway PLS','Multiway CPLS','Multiway CPLS (Yadd)','Location','SouthEast')
pp = get(gcf,'PaperPosition'); set(gcf,'PaperPosition',[0,0,pp(3),pp(4)]); set(gcf,'PaperSize',[pp(3),pp(4)])



%% Three-dimensional X-data
% clear, clc, close all
if(~isempty(dir('C:\Users\kristl')))
    load('C:\Users\kristl\Dropbox\Datasett\Multiway\sugar.mat')
else
    load('C:\Users\ulfin\Dropbox\Dokumenter\Jobb\Arbeider\Datasett\Multiway\sugar.mat')
end
X = reshape(X,DimX);

ncomp = 20;

% Two responses
Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,1:2);
Xval = X(2:2:end,:,:); Yval = y(2:2:end,1:2);
% [B,T,W,Q,R,Wa] = ncpls(ncomp, Xcal, Ycal, [], false);
[B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal);
Tval = GMP(Xval-mean(Xcal,1), R, 2);  % Predicted scores

% Score plot
figure
plot(T(:,1),T(:,2),'o')
hold on, grid on
plot(Tval(:,1),Tval(:,2),'+')
xlabel('Comp. 1'); ylabel('Comp. 2');
legend('Training','Test')

ypred = GMP(Xval-mean(Xcal,1), B, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]); % Predicted responses
RMSEP = zeros(ncomp,size(Yval,2)); R2 = RMSEP; R2cal = RMSEP;
for i=1:ncomp
    RMSEP(i,:) = sqrt(mean((Yval-squeeze(ypred(:,i,:))).^2));
    R2(i,:) = (1 - mean((Yval-squeeze(ypred(:,i,:))).^2)./var(Yval))*100;
    yhat = T(:,1:i)*Q(:,1:i)' + mean(Ycal);             % Fitted values
    R2cal(i,:) = (1 - mean((Ycal-yhat).^2)./var(Ycal))*100;
end

% Two responses & one additional
Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,1:2); YcalAdd = y(1:2:end,3);
Xval = X(2:2:end,:,:); Yval = y(2:2:end,1:2);
[B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal, YcalAdd);
% [B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false);
yhat = T*Q' + mean(Ycal);             % Fitted values
T2 = GMP(Xval-mean(Xcal,1), R, 2);    % Predicted scores
ypred = GMP(Xval-mean(Xcal,1), B, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]); % Predicted responses
RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
for i=1:ncomp
    RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypred(:,i,:))).^2));
    R2c(i,:) = (1 - mean((Yval-squeeze(ypred(:,i,:))).^2)./var(Yval))*100;
end

% Explained variance
figure, hold on
plot(R2cal(:,1), '-k'), plot(R2cal(:,2), '--k')
plot(R2(:,1), '-r'),    plot(R2(:,2), '--r')
plot(R2c(:,1), '-b'),   plot(R2c(:,2), '--b')
ylabel('R^2'), xlabel('#comp')
grid on
legend('Resp. 1, cal.', 'Resp. 2, cal.', 'Resp. 1, val.', 'Resp. 2, val.', 'Resp. 1, Yadd', 'Resp. 2, Yadd','Location','SouthEast')
% legend('Resp. 1, cal.', 'Resp. 2, cal.', 'Resp. 1, val.', 'Resp. 2, val.','Location','SouthEast')


%% Multiway PLS
addpath(genpath('../NPLS'))

[Xfactors,Yfactors,Core,Bm,ypred,ssx,ssy,reg] = npls(Xcal,Ycal,ncomp);
RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
for i=1:ncomp
    [ypreda,Ta,ssXa,Xresa]=npred(Xval,i,Xfactors,Yfactors,Core,Bm);
    RMSEPa(i,:) = sqrt(mean((Yval-(ypreda)).^2));
    R2a(i,:) = (1 - mean((Yval-(ypreda)).^2)./var(Yval))*100;
end

% Explained variance
figure, hold on
plot(R2a(:,1), '-k'), plot(R2a(:,2), '--k')
plot(R2(:,1), '-r'),  plot(R2(:,2), '--r')
plot(R2c(:,1), '-b'), plot(R2c(:,2), '--b')
ylabel('R^2'), xlabel('#comp')
grid on
legend('NPLS, Resp. 1', 'NPLS, Resp. 2', 'NCPLS, Resp. 1', 'NCPLS, Resp. 2', 'NCPLS-add, Resp. 1', 'NCPLS-add, Resp. 2','Location','SouthEast')


%% Direct use of w or outerproduct of processed w
ncomp = 4;
Xcal = X(1:2:end,:,:); Ycal = y(1:2:end,1:2);
Xval = X(2:2:end,:,:); Yval = y(2:2:end,1:2);

% Process w for each component, e.g. [w1,~,w2] = svds(w,1);
[B,T,W,Q,R,Wa] = ncpls(ncomp, Xcal, Ycal);

% Use w without processing
[Bu,Tu,Wu,Qu,Ru,Wau] = ncpls(ncomp, Xcal, Ycal, [], false);

figure
subplot(221)
plot(W{1}); title('Loading weights, internally processed'); axis tight
subplot(222)
plot(Wu{1}); title('Loading weights, post-processed'); axis tight
subplot(223)
plot(W{2}); title('Loading weights, internally processed'); axis tight
subplot(224)
plot(Wu{2}); title('Loading weights, post-processed'); axis tight


%% Milk fats
% clear, clc, close all
load('../Mishra/milk_fats.mat')

ncomp = 50;

% One response
Xcal = X{2}(1:2:end,:,:); Ycal = protein(1:2:end,:);
Xval = X{2}(2:2:end,:,:); Yval = protein(2:2:end,:);
% [B,T,W,Q,R,Wa] = ncpls(ncomp, Xcal, Ycal, [], false);
[B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal);
Tval = GMP(Xval-mean(Xcal,1), R, 2);  % Predicted scores

% Score plot
figure
plot(T(:,1),T(:,2),'o')
hold on, grid on
plot(Tval(:,1),Tval(:,2),'+')
xlabel('Comp. 1'); ylabel('Comp. 2');
legend('Training','Test')

ypred = GMP(Xval-mean(Xcal,1), B, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]); % Predicted responses
RMSEP = zeros(ncomp,size(Yval,2)); R2 = RMSEP; R2cal = RMSEP;
for i=1:ncomp
    RMSEP(i,:) = sqrt(mean((Yval-squeeze(ypred(:,i,:))).^2));
    R2(i,:) = (1 - mean((Yval-squeeze(ypred(:,i,:))).^2)./var(Yval))*100;
    yhat = T(:,1:i)*Q(:,1:i)' + mean(Ycal);             % Fitted values
    R2cal(i,:) = (1 - mean((Ycal-yhat).^2)./var(Ycal))*100;
end

% One response & one additional
Xcal = X{2}(1:2:end,:,:); Ycal = protein(1:2:end,:); YcalAdd = fats(1:2:end,:);
Xval = X{2}(2:2:end,:,:); Yval = protein(2:2:end,:);
[B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal, YcalAdd);
% [B,T,W,Q,R] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false);
yhat = T*Q' + mean(Ycal);             % Fitted values
T2 = GMP(Xval-mean(Xcal,1), R, 2);    % Predicted scores
ypred = GMP(Xval-mean(Xcal,1), B, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]); % Predicted responses
RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
for i=1:ncomp
    RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypred(:,i,:))).^2));
    R2c(i,:) = (1 - mean((Yval-squeeze(ypred(:,i,:))).^2)./var(Yval))*100;
end

% Explained variance
figure, hold on
plot(R2cal(:,1), '-k')
plot(R2(:,1), '-r')
plot(R2c(:,1), '-b')
ylabel('R^2'), xlabel('#comp')
grid on
legend('Cal.', 'Val.', 'Yadd','Location','SouthEast')
% legend('Resp. 1, cal.', 'Resp. 2, cal.', 'Resp. 1, val.', 'Resp. 2, val.','Location','SouthEast')


%% Effect of orthogonalization of w_j and unfolded analysis
% Three-dimensional X-data
% clear, clc, close all
load('../Mishra/milk_fats.mat')

ncomp = 10; %20;

% One response
Xcal = X{2}(1:2:end,:,:); Ycal = protein(1:2:end,:); YcalAdd = fats(1:2:end,:); YcalAdd = [];
Xval = X{2}(2:2:end,:,:); Yval = protein(2:2:end,:);

% Orth w_j, folded analysis
[Ba,Ta,Wa,Qa,Ra,Woa] = ncpls(ncomp, Xcal, Ycal, YcalAdd);
ypreda = GMP(Xval-mean(Xcal,1), Ba, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% Non-orth w_j, folded analysis
[Bb,Tb,Wb,Qb,Rb,Wob] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false);
ypredb = GMP(Xval-mean(Xcal,1), Bb, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% Orth w_j, unfolded analysis
[Bc,Tc,Wc,Qc,Rc,Woc] = ncpls(ncomp, Xcal, Ycal, YcalAdd, true, false);
ypredc = GMP(Xval-mean(Xcal,1), Bc, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);
% Non-orth w_j, unfolded analysis
[Bd,Td,Wd,Qd,Rd,Wod] = ncpls(ncomp, Xcal, Ycal, YcalAdd, false, false);
ypredd = GMP(Xval-mean(Xcal,1), Bd, 2)+reshape(mean(Ycal),[1,1,size(Ycal,2)]);

RMSEPa = zeros(ncomp,size(Yval,2)); R2a = RMSEPa;
RMSEPb = zeros(ncomp,size(Yval,2)); R2b = RMSEPb;
RMSEPc = zeros(ncomp,size(Yval,2)); R2c = RMSEPc;
for i=1:ncomp
    RMSEPa(i,:) = sqrt(mean((Yval-squeeze(ypreda(:,i,:))).^2));
    R2a(i,:) = (1 - mean((Yval-squeeze(ypreda(:,i,:))).^2)./var(Yval))*100;
    RMSEPb(i,:) = sqrt(mean((Yval-squeeze(ypredb(:,i,:))).^2));
    R2b(i,:) = (1 - mean((Yval-squeeze(ypredb(:,i,:))).^2)./var(Yval))*100;
    RMSEPc(i,:) = sqrt(mean((Yval-squeeze(ypredc(:,i,:))).^2));
    R2c(i,:) = (1 - mean((Yval-squeeze(ypredc(:,i,:))).^2)./var(Yval))*100;
end

figure
hold on
plot(R2a(:,1), '-k')
plot(R2b(:,1), '-r')
plot(R2c(:,1), '-b')
ylim([0 100])
grid on
title('Orthogonal w_j')
legend('Orth. w_j','Non-orth w_j','Unfolded', 'Location', 'SouthEast')
xlabel('# comp.')
ylabel('% explained variance')

figure
P = get(gcf,'Position'); set(gcf,'Position', [P(1:2) 560, 145])
subplot(131)
imagesc(Wa{1}'*Wa{1})
title("Orth. w_j : W_1'*W_1")
subplot(132)
imagesc(Wb{1}'*Wb{1})
title("Non-orth. w_j : W_1'*W_1")
subplot(133)
imagesc(Wc{1}'*Wc{1})
title("Unfolded : W_1'*W_1")

for i=1:ncomp
    Wodd1 = GOuter({Wc{1}(:,i),Wc{2}(:,i)});
    Wodd2 = Woc(:,:,i);
    disp(corr(Wodd1(:),Wodd2(:)))
    disp(norm(Wodd2-Wodd1, 'fro'))
end


%% Direct use of w or outerproduct of processed w
ncomp = 4;
Xcal = X{2}(1:2:end,:,:); Ycal = protein(1:2:end,:);
Xval = X{2}(2:2:end,:,:); Yval = protein(2:2:end,:);

% Process w for each component, e.g. [w1,~,w2] = svds(w,1);
[B,T,W,Q,R,Wa] = ncpls(ncomp, Xcal, Ycal);

% Use w without processing
[Bu,Tu,Wu,Qu,Ru,Wau] = ncpls(ncomp, Xcal, Ycal, [], false);

figure
subplot(221)
plot(W{1}); title('Loading weights, internally processed'); axis tight
subplot(222)
plot(Wu{1}); title('Loading weights, post-processed'); axis tight
subplot(223)
plot(W{2}); title('Loading weights, internally processed'); axis tight
subplot(224)
plot(Wu{2}); title('Loading weights, post-processed'); axis tight

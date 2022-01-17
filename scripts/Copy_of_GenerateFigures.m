clearvars;

outdir = '../results/IgorFiles/';
datadir = '../data/';

outputIgor = true;

%%
% AnyAMD_ActiveWetAMD figure

% Load the CSV file.
% fid = load([datadir, 'AnyAMD_ActiveWetAMD.csv']);

T = readtable([datadir, 'AnyAMD_ActiveWetAMD.csv']);

% Loop through the table and grab the columns.
labels = T.Var1; % Predictor
% Coefficient
coeff = T.Var2;
% CI lower
ciLower = T.Var3;
% CI upper
ciUpper = T.Var4;

x = 1:(length(coeff)-1);

redPts = (ciLower > 1 | ciUpper < 1);

colors = [
    0,0,0
    0.8,0,0 
    ];

condLabels = {
    'NonSig'
    'Sig'
    };

%
figure(30); clf;

subplot(221)
hold on
for j = 1 : 2
    pts = (redPts == j-1);
    M = nan(size(coeff));
    L = nan(size(coeff));
    U = nan(size(coeff));
    M(pts) = coeff(pts);
    L(pts) = -(ciLower(pts)-coeff(pts));
    U(pts) = ciUpper(pts)-coeff(pts);
    
    line(M(2:end), x, 'DisplayName', [condLabels{j},'Mean'], 'Color', colors(j,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    line(L(2:end), x, 'DisplayName', [condLabels{j},'Lower'], 'Color', colors(j,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    line(U(2:end), x, 'DisplayName', [condLabels{j},'Upper'], 'Color', colors(j,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
end
line(ones(1,2), x([1,end]), 'DisplayName', 'line', 'Color', 'k', 'LineStyle', '--', 'Marker', 'none', 'LineWidth', 2);
hold off;
xlabel('predictor'); ylabel('value');
if outputIgor
    fname = 'AnyAMD_ActiveWetAMD';
    makeAxisStruct(gca, fname, 'basedir', outdir);
    hdf5write([outdir,fname,'.h5'], [fname, '/labels_Y'], x, 'WriteMode', 'append');
    hdf5write([outdir,fname,'.h5'], [fname, '/labels_X'], labels(2:end), 'WriteMode', 'append');
end




return
figure(30); clf;
hold on
for j = 2 : length(coeff)
%     line(x(j-1),coeff(j), 'DisplayName', ['x',num2str(x(j-1)),'Mean'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
%     line(x(j-1),ciLower(j)-coeff(j), 'DisplayName', ['x',num2str(x(j-1)),'L'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
%     line(x(j-1),ciUpper(j)-coeff(j), 'DisplayName', ['x',num2str(x(j-1)),'U'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    line(coeff(j), x(j-1), 'DisplayName', ['x',num2str(x(j-1)),'Mean'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    line(-(ciLower(j)-coeff(j)), x(j-1), 'DisplayName', ['x',num2str(x(j-1)),'L'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    line(ciUpper(j)-coeff(j), x(j-1), 'DisplayName', ['x',num2str(x(j-1)),'U'], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
    
%     line(coeff(j), x(j-1), 'DisplayName', ['coeff',num2str(x(j-1))], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
%     line(ciLower(j), x(j-1), 'DisplayName', ['lower',num2str(x(j-1))], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
%     line(ciUpper(j), x(j-1), 'DisplayName', ['upper',num2str(x(j-1))], 'Color', colors(redPts(j)+1,:), 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 2);
end
line(ones(1,2), x([1,end]), 'DisplayName', 'line', 'Color', colors(redPts(j)+1,:), 'LineStyle', '--', 'Marker', 'none', 'LineWidth', 2);
hold off;
xlabel('predictor'); ylabel('value');
if outputIgor
    fname = 'AnyAMD_ActiveWetAMD';
    makeAxisStruct(gca, fname, 'basedir', outdir);
    hdf5write(fileName, strcat(dataRoot, '/', ops.prefix, fNames{i}), s.(fNames{i}), 'WriteMode', cWrMode);
end

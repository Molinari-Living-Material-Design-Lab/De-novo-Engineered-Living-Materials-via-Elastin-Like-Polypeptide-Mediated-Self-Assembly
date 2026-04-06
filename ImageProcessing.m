clear; clc; close all;

alanineDir = "C:\Users\User\MATLAB Drive\ImageAnalysis_Rotation\Alanine10_Images";
glutDir    = "C:\Users\User\MATLAB Drive\ImageAnalysis_Rotation\GlutamicAcid10_Images";
lysineDir  = "C:\Users\User\MATLAB Drive\ImageAnalysis_Rotation\Lysine10_Images";

alanineFiles = {
    'Alanine10_1_10292025.JPG'
    'Alanine10_2_10292025.JPG'
    'Alanine10_3_10292025.JPG'
    'Alanine10_4_10292025.JPG'
    'Alanine10_5_10292025.JPG'
    'Alanine10_7_10292025.JPG'
    'Alanine10_8_10292025.JPG'
    'Alanine10_9_10292025.JPG'
    'Alanine10_10_10292025.JPG'
    'Alanine10_12_10292025.JPG'
    'Alanine10_16_10292025.JPG'
};

glutFiles = {
    'GlutamicAcid_1_10302025_New.JPG'
    'GlutamicAcid_2_10302025_New.JPG'
    'GlutamicAcid_3_10302025_New.JPG'
    'GlutamicAcid_4_10302025_New.JPG'
    'GlutamicAcid_6_10302025_New.JPG'
    'GlutamicAcid_8_10302025_Old.JPG'
    'GlutamicAcid_9_10302025_Old.JPG'
    'GlutamicAcid_10_10302025_Old.JPG'
    'GlutamicAcid_11_10302025_Old.JPG'
    'GlutamicAcid_12_10302025_Old.JPG'
};

lysineFiles = {
    'Lysine10_1_10282025.JPG'
    'Lysine10_2_10282025.JPG'
    'Lysine10_4_10282025.JPG'
    'Lysine10_5_10282025.JPG'
    'Lysine10_6_10282025.JPG'
    'Lysine10_7_10282025.JPG'
    'Lysine10_8_10282025.JPG'
    'Lysine10_9_10282025.JPG'
    'Lysine10_10_10282025.JPG'
    'Lysine10_11_10282025.JPG'
    'Lysine10_12_10282025.JPG'
};

alanineFull = fullfile(alanineDir, string(alanineFiles));
glutFull    = fullfile(glutDir,    string(glutFiles));
lysineFull  = fullfile(lysineDir,  string(lysineFiles));
imageFiles  = [alanineFull; glutFull; lysineFull];
numFiles    = numel(imageFiles);

if numFiles == 0
    error('No files listed.');
end

roiShrink = 0.95;
innerMarginPx = 6;

minPlateArea = 10000;
minFilamentArea = 15;
closeRad = 2;
bridgeOn = true;
dilateRad = 1;

otsuScale = 0.65;
useAdaptiveFallback = true;
adaptSensitivity = 0.48;
adaptWindowFrac = 0.08;

stdGate = 0.010;

ROI_DIAMETER_MM = 67;
ROI_AREA_MM2 = pi*(ROI_DIAMETER_MM/2)^2;

outFolder = fullfile(pwd, "Results_HSV");
if ~exist(outFolder, "dir"), mkdir(outFolder); end

ImageName        = strings(0,1);
SetName          = strings(0,1);
Status           = strings(0,1);
ROIDiameter_px   = zeros(0,1);
ROIArea_px2      = zeros(0,1);
FilamentCov_pct  = zeros(0,1);
FilamentArea_px2 = zeros(0,1);
FilamentArea_mm2 = zeros(0,1);

for k = 1:numFiles
    fullFileName = imageFiles(k);
    [~, baseFileName, ext] = fileparts(fullFileName);
    baseFileNameWithExt = baseFileName + ext;

    if startsWith(baseFileName, "Alanine10_")
        setNameThis = "Alanine";
    elseif startsWith(baseFileName, "GlutamicAcid_")
        setNameThis = "GlutamicAcid";
    elseif startsWith(baseFileName, "Lysine10_")
        setNameThis = "Lysine";
    else
        setNameThis = "Unknown";
    end

    if ~isfile(fullFileName)
        fprintf('Missing [%d/%d]: %s\n', k, numFiles, fullFileName);
        fprintf('------------------------------------------\n');
        continue;
    end

    fprintf('Processing [%d/%d]: %s\n', k, numFiles, fullFileName);

    img = imread(fullFileName);
    imgD = im2double(img);

    hsv = rgb2hsv(img);
    h = hsv(:,:,1);
    s = hsv(:,:,2);
    v = hsv(:,:,3);

    yellowMask = (h > 0.08 & h < 0.25) & (s > 0.12) & (v > 0.10);
    yellowMask = imfill(yellowMask, 'holes');
    yellowMask = bwareaopen(yellowMask, minPlateArea);

    props = regionprops(yellowMask, 'Centroid', 'EquivDiameter', 'Area', 'BoundingBox');

    if isempty(props)
        fprintf('No plate detected in %s.\n', baseFileNameWithExt);
        fprintf('------------------------------------------\n');
        continue;
    end

    [~, idx] = max([props.Area]);
    center = props(idx).Centroid;
    radius = (props(idx).EquivDiameter/2) * roiShrink;

    [cols, rows] = meshgrid(1:size(img,2), 1:size(img,1));
    roiMask = (cols - center(1)).^2 + (rows - center(2)).^2 <= radius^2;

    inner = imerode(roiMask, strel('disk', innerMarginPx));
    if nnz(inner) < 1000
        inner = roiMask;
    end

    bbox = props(idx).BoundingBox;
    buf = 50;
    yRange = max(1, floor(bbox(2)-buf)) : min(size(img,1), ceil(bbox(2)+bbox(4)+buf));
    xRange = max(1, floor(bbox(1)-buf)) : min(size(img,2), ceil(bbox(1)+bbox(3)+buf));

    lab = rgb2lab(imgD);
    L = lab(:,:,1);
    Ln = mat2gray(L);
    darkness = 1 - Ln;

    target = 0.55*s + 0.45*darkness;
    target(~inner) = 0;

    seRad = max(8, round(0.06 * radius));
    se = strel('disk', seRad);

    background = imopen(target, se);
    itop = imsubtract(target, background);
    itop = mat2gray(itop);
    itop(~inner) = 0;

    roiValues = itop(inner);
    stdDev = std(double(roiValues));

    if stdDev < stdGate || isempty(roiValues)
        filamentMask = false(size(img,1), size(img,2));
        statusStr = "Clean Plate (low contrast)";
    else
        level = graythresh(roiValues);
        bwOtsu = imbinarize(itop, level * otsuScale);
        bwOtsu = bwOtsu & inner;

        if useAdaptiveFallback
            frac = nnz(bwOtsu) / max(nnz(inner),1);
            if frac < 0.0005
                win = max(15, round(adaptWindowFrac * (2*radius)));
                if mod(win,2)==0, win = win+1; end
                T = adaptthresh(itop, adaptSensitivity, ...
                    'NeighborhoodSize', [win win], 'ForegroundPolarity','bright');
                bwA = imbinarize(itop, T) & inner;
                filamentMask = bwA;
                statusStr = "Filaments Detected (adaptive fallback)";
            else
                filamentMask = bwOtsu;
                statusStr = "Filaments Detected (Otsu)";
            end
        else
            filamentMask = bwOtsu;
            statusStr = "Filaments Detected (Otsu)";
        end
    end

    filamentMask = bwareaopen(filamentMask, minFilamentArea);
    filamentMask = imclose(filamentMask, strel('disk', closeRad));

    if bridgeOn
        filamentMask = bwmorph(filamentMask, 'bridge');
    end
    if dilateRad > 0
        filamentMask = imdilate(filamentMask, strel('disk', dilateRad));
    end

    filamentMask = filamentMask & roiMask;

    roiAreaPX2 = nnz(roiMask);
    filamentAreaPX2 = nnz(filamentMask);
    coveragePCT = (filamentAreaPX2 / max(roiAreaPX2,1)) * 100;
    roiDiamPX = 2*sqrt(roiAreaPX2/pi);

    mm2_per_px2 = ROI_AREA_MM2 / max(roiAreaPX2,1);
    filamentAreaMM2 = filamentAreaPX2 * mm2_per_px2;

    fprintf('Status: %s\n', statusStr);
    fprintf('ROI diameter (px): %.2f\n', roiDiamPX);
    fprintf('ROI area (px^2): %d\n', roiAreaPX2);
    fprintf('Filament coverage (%%): %.4f\n', coveragePCT);
    fprintf('Filament area (px^2): %d\n', filamentAreaPX2);
    fprintf('Filament area (mm^2): %.4f\n', filamentAreaMM2);

    ImageName(end+1,1)        = baseFileNameWithExt;
    SetName(end+1,1)          = setNameThis;
    Status(end+1,1)           = statusStr;
    ROIDiameter_px(end+1,1)   = roiDiamPX;
    ROIArea_px2(end+1,1)      = roiAreaPX2;
    FilamentCov_pct(end+1,1)  = coveragePCT;
    FilamentArea_px2(end+1,1) = filamentAreaPX2;
    FilamentArea_mm2(end+1,1) = filamentAreaMM2;

    figure(k);
    set(gcf, 'Name', sprintf('File is %s', baseFileNameWithExt), ...
        'Units', 'Normalized', 'Position', [0.1 0.1 0.7 0.4]);

    imgZoom = img(yRange, xRange, :);
    filamentMaskZoom = filamentMask(yRange, xRange);

    subplot(1,2,1);
    imshow(img); hold on;
    theta = linspace(0, 2*pi, 200);
    plot(center(1) + radius*cos(theta), center(2) + radius*sin(theta), 'g', 'LineWidth', 2);
    title(sprintf('ROI: %s', baseFileNameWithExt), 'Interpreter', 'none');

    subplot(1,2,2);
    imshow(imgZoom); hold on;
    redOverlay = cat(3, ones(size(imgZoom,1), size(imgZoom,2)), ...
                        zeros(size(imgZoom,1), size(imgZoom,2)), ...
                        zeros(size(imgZoom,1), size(imgZoom,2)));
    hRed = imshow(redOverlay);
    set(hRed, 'AlphaData', filamentMaskZoom * 0.5);
    title(sprintf('Coverage: %.4f%%', coveragePCT));

    drawnow;
    fprintf('------------------------------------------\n');
end

T = table(ImageName, SetName, Status, ROIDiameter_px, ROIArea_px2, FilamentCov_pct, FilamentArea_px2, FilamentArea_mm2);

csvPath = fullfile(outFolder, "FlaskFilaments_Results.csv");
writetable(T, csvPath);

xlsxPath = fullfile(outFolder, "FlaskFilaments_Results.xlsx");
try
    writetable(T, xlsxPath, "FileType", "spreadsheet");
catch
    warning('Could not write XLSX (maybe missing Excel support). CSV saved anyway.');
end

fprintf('Saved CSV:  %s\n', csvPath);
fprintf('Saved XLSX: %s\n', xlsxPath);
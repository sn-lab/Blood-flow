clear
% to calculate average diameter of a vessel inside ROI
% 8/13/2015 ~mh973

%Input
threshPrctile=70;
selectROI=true;
selectImages=true;

%Script
if selectImages
    pname = [uigetdir '/'];
    
end
files = dir([pname '*.tif*']);
txtFile=fopen([pname 'average diameter.txt'],'a');
for i=1:length(files)
    im=imread([pname files(i).name]);
    bw=im2bw(im,single(prctile(im(:),threshPrctile))/ ...
        single(intmax(class(im))));
    if selectROI
        [maskBW,xi,yi]=roipoly(255*bw);
        save(strcat(num2str(i),'.mat'),'maskBW','xi','yi');
    else
        load(strcat(num2str(i),'.mat'),'maskBW','xi','yi');
    end
    bw(~maskBW)=false;
    avgDiameter= sum(bw(:))/pdist([mean(xi(1:2)),mean(yi(1:2)); ...
        mean(xi(3:4)),mean(yi(3:4))])
    fprintf(txtFile,'%s , %f , %d \n',files(i).name,avgDiameter,...
        threshPrctile);
end
    
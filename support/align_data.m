function alignedData = align_data(unalignedData,unalignedTime,alignedTime)
%FUNCTION alignedData = align_data(unalignedData,unalignedTime,alignedTime)
%
% Takes a (1xN) data vector and associated (1xN) time vector, and aligns it 
% to a new (regularly spaced) (1xM) time vector.
%
%INPUTS
% unalignedData: data vector
% unalignedTime: vector of timepoints corresponding to the data vector
% alignedTime: vector of timepoints to align the data to

%check that inputs are valid
assert(size(unalignedTime,1)==1 & size(alignedTime,1)==1 & size(unalignedData,1)==1,...
    'all inputs must be row vectors');
assert(ismatrix(unalignedTime) & ismatrix(alignedTime) & ismatrix(unalignedData),...
    'all inputs must be row vectors');
assert(length(unalignedTime)==length(unalignedData),['unalignedTime and ',...
    'unalignedData vectors must be the same size']);

%calculate the desired sampling rate
timeBetweenAlignedSamples = diff(alignedTime);
alignedSampleRate = median(timeBetweenAlignedSamples);

%check that the desired timepoints are regularly spaced
assert(all((timeBetweenAlignedSamples-alignedSampleRate)<0.001),['alignedTime must '...
    'have a regularly-spaced sampling rate (e.g. [0 .1 .2 .3], not [0 .09 .21 .3]']);

%get indices of alignedTime
alignedIndices = round(0.25+(alignedTime/alignedSampleRate));
alignedColSub = alignedIndices-min(alignedIndices)+1;

%get indices of unalignedTime
unalignedIndices = round(unalignedTime/alignedSampleRate);
unalignedColSub = unalignedIndices-min(alignedIndices)+1;

%only include unalignedTime timepoints that are within the bounds of alignedTime
outOfBoundsInds = unalignedColSub<=0 | unalignedColSub>max(alignedColSub);
unalignedColSub = unalignedColSub(~outOfBoundsInds);
unalignedData = unalignedData(~outOfBoundsInds);
numDatapoints = length(unalignedColSub);

%distribute unalignedData into the appropriate alignedTime time bins (ie columns)
% alignedData = nan(numDatapoints,max(alignedColSub));
% alignedData(sub2ind(size(alignedData),1:numDatapoints,unalignedColSub(s:e))) = unalignedData(s:e);
% 
% %collapse to a row vector (if multiple unaligned timepoints are in the same aligned time bin, they are averaged)
% alignedData = mean(alignedData,'omitnan');

%same as what's commented out above, but this way uses significantly less memory
alignedData = sparse(1:numDatapoints,unalignedColSub,unalignedData,numDatapoints,max(alignedColSub),numDatapoints);
alignedDataOnes = sparse(1:numDatapoints,unalignedColSub,ones(1,numDatapoints),numDatapoints,max(alignedColSub),numDatapoints);

%collapse to a row vector (if multiple unaligned timepoints are in the same aligned time bin, they are averaged)
alignedData = sum(alignedData,'omitnan')./sum(alignedDataOnes,'omitnan');
alignedData = full(alignedData); %de-sparse the data;

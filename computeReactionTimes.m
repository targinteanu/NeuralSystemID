function [rtTable, medianRT, meanRT] = computeReactionTimes(events)
% computeReactionTimes: Compute reaction times from events table
% Inputs:
%   events - table from parseNEVSerialCodes 
% Outputs:
%   rtTable  - table with columns: StartTime, EndTime, ReactionTime_s, EndCode
%   medianRT - median reaction time (seconds)
%   meanRT   - mean reaction time (seconds)
% Notes:
%   - Reaction time is defined as the time difference between decoding start (code 7)
%     and the next decoding end (code 8 or 9) event.
%   - Only pairs where start precedes end are included.

startIdx = find(events.Code == 7);
endIdx = find(events.Code == 8 | events.Code == 9);
% OR endIdx = find(events.Code == 9); gets dif values

rtList = [];
startTimes = [];
endTimes = [];
endCodes = [];

for i = 1:numel(startIdx)
    sIdx = startIdx(i);
    sTime = events.Time_s(sIdx);
    % Find the first end event after this start
    eIdx = endIdx(endIdx > sIdx);
    if ~isempty(eIdx)
        eIdx = eIdx(1);
        eTime = events.Time_s(eIdx);
        rt = eTime - sTime;
        rtList(end+1,1) = rt;
        startTimes(end+1,1) = sTime;
        endTimes(end+1,1) = eTime;
        endCodes(end+1,1) = events.Code(eIdx);
    end
end


if isempty(rtList)
    medianRT = NaN;
    meanRT = NaN;
    medianCol = nan(size(rtList));
    meanCol = nan(size(rtList));
else
    medianRT = median(rtList);
    meanRT = mean(rtList);
    medianCol = repmat(medianRT, size(rtList));
    meanCol = repmat(meanRT, size(rtList));
end

rtTable = table(startTimes, endTimes, rtList, endCodes, medianCol, meanCol, ...
    'VariableNames', {'StartTime_s','EndTime_s','ReactionTime_s','EndCode','MedianRT','MeanRT'});

end

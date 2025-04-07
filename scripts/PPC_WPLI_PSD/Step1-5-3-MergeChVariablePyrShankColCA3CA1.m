%% Combine the best CA1 channels and CA3 channels

%% CA1 Channels
clear
LoadData='Y:\singer\LuZhang\Project2ChronicF\Results\step1\';

load([LoadData 'ThetaAlign\ShankColAlignPyrOriCh.mat'],'ThetaList','ProbeXY','ChsortII','ChsortSh');
CA1List=ThetaList;

%% CA3 channels
CA1List=rmfield(ThetaList,{'CA3ThChIndex','CA3PyrChIndex'});
load([LoadData 'ThetaAlign\CA3ShankColAlignPyrOriCh.mat'],'ThetaList');
CA3List=ThetaList;

for i=1:length(CA1List)
    CA1List(i).CA3ThChIndex=CA3List(i).CA3ThChIndex;
    CA1List(i).CA3PyrChIndex=CA3List(i).CA3PyrChIndex;
end

ThetaList=CA1List;
save([LoadData 'ThetaAlign\ShankColAlignCA1PyrOriCA3PyrCh.mat'],'ThetaList','ProbeXY','ChsortII','ChsortSh');

%%
clear;
close all;

%%
dataWastewater=readtimetable('wastewater_santa_clara_completed.csv');

Gilroy=dataWastewater(strcmp(dataWastewater.city,'Gilroy'),:);
PaloAlto=dataWastewater(strcmp(dataWastewater.city,'Palo Alto'),:);
SanJose=dataWastewater(strcmp(dataWastewater.city,'San Jose'),:);
Sunnyvale=dataWastewater(strcmp(dataWastewater.city,'Sunnyvale'),:);

%%
TR=timerange('1-Dec-2020','1-May-2021');
Gilroy=Gilroy(TR,:);
PaloAlto=PaloAlto(TR,:);
SanJose=SanJose(TR,:);
Sunnyvale=Sunnyvale(TR,:);

%%
viralRNAGilroy=Gilroy.N_Gene_gc_g_dry_weight;
viralRNAPaloAlto=PaloAlto.N_Gene_gc_g_dry_weight;
viralRNASanJose=SanJose.N_Gene_gc_g_dry_weight;
viralRNASunnyvale=Sunnyvale.N_Gene_gc_g_dry_weight;

%%
plot(viralRNAGilroy);
hold on;
plot(viralRNAPaloAlto);
plot(viralRNASanJose);
plot(viralRNASunnyvale);
hold off;
xlabel('Days')
title('Viral RNA (copies/grams of dry sludge): 12/01/2020 to 05/01/2021')
legend('Gilroy','Palo Alto','San Jose', 'Sunnyvale')
xlim([0 100])

%%
dataTests=readtimetable('case_counts.csv');
TotalCases=table2array(dataTests(TR,{'Total_cases'}));

%%
subplot(4,1,1)
plot(log10(viralRNAGilroy));
hold on;
plot(log10(diff(TotalCases)));
hold off;
xlabel('Days')
legend('Gilroy Viral RNA','Positive Results')

subplot(4,1,2)
plot(log10(viralRNAPaloAlto));
hold on;
plot(log10(diff(TotalCases)));
hold off;
xlabel('Days')
legend('Palo Alto Viral RNA','Positive Results')

subplot(4,1,3)
plot(log10(viralRNASanJose));
hold on;
plot(log10(diff(TotalCases)));
hold off;
xlabel('Days')
legend('San Jose Viral RNA','Positive Results')

subplot(4,1,4)
plot(log10(viralRNASunnyvale));
hold on;
plot(log(diff(TotalCases)));
hold off;
xlabel('Days')
legend('Sunnyvale Viral RNA','Positive Results')

%%
% save('viralRNAGilroy.mat','Gilroy')
% save('viralRNAPaloAlto.mat','PaloAlto')
% save('viralRNASanJose.mat','SanJose')
% save('viralRNASunnyvale.mat','Sunnyvale')
% save('TotalCases.mat','TotalCases')

%%
% [r,lag]=xcorr(TotalCases,viralRNAGilroy);
% [m,idx]=max(r); 
% plot(lag,r)
% xline(lag(idx),'--')
% xlabel('Lag')
% ylabel('Correlation')
% title(['Lag Time: ' num2str(idx-150)])









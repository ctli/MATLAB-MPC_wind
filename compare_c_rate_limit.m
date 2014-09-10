clear all
close all
clc

% alpha = 0.79; % weighting coefficient in the persistent algorithm

% no c rate limits
a = [-917339.7445	-1421308.683	-1651443.553	-1718233.854	-1827200.07	-1868194.137
-981243.323	-1545000.867	-1799601.573	-1873638.654	-2010304.624	-2086994.702
];

% w/ c rate limits
b = [-394117.1509	-581590.3575	-757142.4981	-922469.2864	-1427949.85	-1866653.667
-981243.323	-1545000.867	-1799601.573	-1873638.654	-2010304.624	-2086994.702
];

%%
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);
fcst = data.data(:,4);

scale = 800/max(fcst);
obs = obs*scale;
fcst = fcst*scale;

obs = obs(1:8760)';
fcst = fcst(1:8760)';

wind_nameplate = 800; % [MW]

%%
c1 = 1;
c2 = 1.03;
c3 = 1;
full_potential = -c1*sum(obs);


%%
load ConvReserve
total_cost_ConvReserve = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched);
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve = sum(temp);

%%
figure(1); clf; hold on;
line([0.5 7.5], [1 1]*0.2, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.4, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.6, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.8, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*1, 'color', [0.9 0.9 0.9]);

h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'edge', 'none'); hold on;
h2 = bar(2:6+1, a'/full_potential, 1, 'group', 'edge', 'none');
set(h2(1), 'facec', [1 0 0]);
set(h2(2), 'facec', [255 110 110]/255);

h3 = bar(2:6+1, b'/full_potential, 0.5, 'group', 'edge', 'none');

legend([h2, h3], '\eta=0.75, no c lmt', '\eta=0.95, no c lmt', '\eta=0.75, w/ c lmt', '\eta=0.95, w/ c lmt', 2);

xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 8, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});


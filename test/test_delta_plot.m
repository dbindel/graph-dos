% Illustrate plot of a filter polynomial centered at 0.5
%
c = moments_delta(0.5,200);
c(1) = c(1)/2;
cf = filter_jackson(c);

figure;
hold on;
subplot(2,1,1); plot_chebp(c); title('Unfiltered')
subplot(2,1,2); plot_chebp(cf); title('Filtered')
hold off;

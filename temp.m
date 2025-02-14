
skew0 = abs(v0(:,2));
s01 = v0(:,3);
s02 = v0(:,4);

skew1 = abs(v1(:,2));
s11 = v1(:,3);
s12 = v1(:,4);

figure(1);
subplot(2,2,1)
yyaxis left;
plot(log(skew0),s01,'bo');
subplot(2,2,3)
yyaxis right;
plot(log(skew0),s02,'ro'); ylim([min(s02) max(s02)]);

subplot(2,2,2)
yyaxis left;
plot(log(skew1),s11,'bo'); 
subplot(2,2,4)
yyaxis right;
plot(log(skew1),s12,'ro'); ylim([min(s12) max(s12)]);

%figure(2);
%yyaxis left;
%plot(log(skew0),s02,'bo', 'MarkerFaceColor','b');
% yyaxis right;
% plot(log(skew1),s12,'ro', 'MarkerFaceColor','r');

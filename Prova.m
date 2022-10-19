clc
load S818_10.txt
load S831_10.txt
load S832_10.txt
load S827_10.txt
load S828_10.txt
load S830_10.txt

plot(S832_10(:,1),S832_10(:,2))
hold on
plot(S828_10(:,1),S828_10(:,2))
hold on
plot(S831_10(:,1),S831_10(:,2))
plot(S830_10(:,1),S830_10(:,2))
plot(S827_10(:,1),S827_10(:,2))

plot(S818_10(:,1),S818_10(:,2))
legend('S832','S828','S831','S818','S830','S827')
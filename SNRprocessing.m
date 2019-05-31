S=[max(1-A0(end/2:end)),max(1-A5(end/2:end)),max(1-A10(end/2:end)),max(1-A15(end/2:end)),max(1-A20(end/2:end)),max(1-A25(end/2:end)),max(1-A30(end/2:end)),max(1-A35(end/2:end)),max(1-A40(end/2:end)),max(1-A45(end/2:end)),max(1-A50(end/2:end)),max(1-A55(end/2:end)),max(1-A60(end/2:end)),max(1-A65(end/2:end)),max(1-A70(end/2:end)),max(1-A75(end/2:end)),max(1-A80(end/2:end)),max(1-A85(end/2:end)),max(1-A90(end/2:end)),max(1-A95(end/2:end)),max(1-A100(end/2:end))];
%M=[1-A0(end),1-A5(end),1-A10(end),1-A15(end),1-A20(end),1-A25(end),1-A30(end),1-A35(end),1-A40(end),1-A45(end),1-A50(end),1-A55(end),1-A60(end),1-A65(end),1-A70(end),1-A75(end),1-A80(end),1-A85(end),1-A90(end),1-A95(end),1-A100(end)];

l=0:5:100;
plot(l,S./M);
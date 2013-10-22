ox_gs_mchg = zeros(1,length(ox_zm));
for i = 1:length(oxazine.aorbAtom)
    % disp([num2str(i),', ',num2str(oxazine.aorbAtom(i)),', ',num2str(ox_gsdm(i,i))]);
ox_gs_mchg(oxazine.aorbAtom(i)) = ox_gs_mchg(oxazine.aorbAtom(i)) - ox_gsdm(i,i);
end

ox_ex_mchg = zeros(1,length(ox_zm));
for i = 1:length(oxazine.aorbAtom)
    % disp([num2str(i),', ',num2str(oxazine.aorbAtom(i)),', ',num2str(ox_gsdm(i,i))]);
ox_ex_mchg(oxazine.aorbAtom(i)) = ox_ex_mchg(oxazine.aorbAtom(i)) - ox_exdm(i,i);
end

rf_gs_mchg = zeros(1,length(rf_zm));
for i = 1:length(resorufin.aorbAtom)
    % disp([num2str(i),', ',num2str(oxazine.aorbAtom(i)),', ',num2str(ox_gsdm(i,i))]);
rf_gs_mchg(resorufin.aorbAtom(i)) = rf_gs_mchg(resorufin.aorbAtom(i)) - rf_gsdm(i,i);
end

rf_ex_mchg = zeros(1,length(rf_zm));
for i = 1:length(resorufin.aorbAtom)
    % disp([num2str(i),', ',num2str(oxazine.aorbAtom(i)),', ',num2str(ox_gsdm(i,i))]);
rf_ex_mchg(resorufin.aorbAtom(i)) = rf_ex_mchg(resorufin.aorbAtom(i)) - rf_exdm(i,i);
end

for i = 1:length(ox_zm)
switch lower(ox_zm{i,2})
case 'c'
ox_gs_mchg(i) = ox_gs_mchg(i) + 4;
ox_ex_mchg(i) = ox_ex_mchg(i) + 4;
case 'n'
ox_gs_mchg(i) = ox_gs_mchg(i) + 5;
ox_ex_mchg(i) = ox_ex_mchg(i) + 5;
case 'h'
ox_gs_mchg(i) = ox_gs_mchg(i) + 1;
ox_ex_mchg(i) = ox_ex_mchg(i) + 1;
case 'o'
ox_gs_mchg(i) = ox_gs_mchg(i) + 6;
ox_ex_mchg(i) = ox_ex_mchg(i) + 6;
end
end

for i = 1:length(rf_zm)
switch lower(rf_zm{i,2})
case 'c'
rf_gs_mchg(i) = rf_gs_mchg(i) + 4;
rf_ex_mchg(i) = rf_ex_mchg(i) + 4;
case 'n'
rf_gs_mchg(i) = rf_gs_mchg(i) + 5;
rf_ex_mchg(i) = rf_ex_mchg(i) + 5;
case 'h'
rf_gs_mchg(i) = rf_gs_mchg(i) + 1;
rf_ex_mchg(i) = rf_ex_mchg(i) + 1;
case 'o'
rf_gs_mchg(i) = rf_gs_mchg(i) + 6;
rf_ex_mchg(i) = rf_ex_mchg(i) + 6;
end
end
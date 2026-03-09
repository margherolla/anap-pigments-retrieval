%% extract pigments from gaussians peaks

function [pigm_der] = derive_pigm_gaus(agaus, stations, A, B)


for i = 1:length(stations)
    agaus_st = agaus(i, :);
    for j = 1:size(agaus, 2)
    pigm_der (i, j) = ((agaus_st(j)./10^(A(j))).^(1/B(j)));
    end
end

end
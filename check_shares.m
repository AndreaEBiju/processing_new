S = load('real_layer12_results.mat');
pc = S.resultsPerChannel.perChannelDetail;
tf = pc.RVN.thetaFull;
uM = [10 50 100 10 50 100 10 50 100];
uE = [10 10 10 100 100 100 1000 1000 1000];
[p_cap, p_1, p_2] = model_layer12_equations('term_shares', uM, uE, tf);
for i = 1:numel(uM)
    fprintf('M%d E%d: p_cap=%.6f p_1=%.6e p_2=%.6e\n', uM(i), uE(i), p_cap(i), p_1(i), p_2(i));
end

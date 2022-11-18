function[acertos_training_pct, alphabar] = MinimosQuadrados_Training(treinamento)
    ytraining = treinamento(:, 11);
    Xtraining = [ones(ytraining), treinamento(:, 1:10)];
    alphabar = Gaussian_Elimination_4(Xtraining' * Xtraining, Xtraining' * ytraining);
    prevision = Xtraining * alphabar;
    conf_prevision = prevision .* ytraining;
    acertos_training_pct = (sum(conf_prevision >= 0) / size(Xtraining, "r"));
endfunction

function[acertos_testing_pct] = MinimosQuadrados_Testing(testagem, alphabar)
    ytesting = testagem(:, 11);
    Xtesting = [ones(ytesting), testagem(:, 1:10)];
    prevision = Xtesting * alphabar;
    conf_prevision = prevision .* ytesting;
    acertos_testing_pct = (sum(conf_prevision >= 0) / size(Xtesting, "r"));
endfunction


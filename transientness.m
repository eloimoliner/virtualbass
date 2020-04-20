function Y = transientness(X,nMedianH,nMedianV)

%% Median filtering
%Vertical
X_v_median = medfilt1(abs(X),nMedianV,[],1);
%Horizontal
X_h_median = medfilt1(abs(X),nMedianH,[],2);

% Compute transientness
Y = (X_v_median.^2) ./ (X_v_median.^2 + X_h_median.^2);
Y(isnan(Y)) = 0;

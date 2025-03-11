% Function for linear fitting with a zero y-intercept.
function slope = linearfit(x, y)
    % Ensure x and y are column vectors
    if isrow(x), x = x' ; end
    if isrow(y), y = y' ; end
    
    % Compute the slope using least squares (Y = slope*X)
    slope = (x' * y) / (x' * x); % The best-fit slope minimizing the sum of squared errors.    
end

function Sigma_thresholded = soft_thresholding1(Sigma,lambda)
    % Soft-thresholding function Y_lambda(sigma)

    % Initialize the thresholded singular values
    Sigma_thresholded = zeros(size(Sigma));
    
    for i = 1:length(Sigma)
        sigma_i = Sigma(i);
         
        % Calculate the threshold based on the given lambda
        threshold = (lambda / 8) * (abs(sigma_i) /3)^(-2/3);
        
        if abs(sigma_i) > (nthroot(54, 3) / 4) * (lambda^(2/3))
            % Apply the soft-thresholding operation
            Sigma_thresholded(i) = (2/3) * sigma_i * (1 + cos(2 * pi / 3 - 2 / 3 * acos(threshold)));
        else
            % Set the value to zero if it does not meet the condition
            Sigma_thresholded(i) = 0;
        end
    end
end
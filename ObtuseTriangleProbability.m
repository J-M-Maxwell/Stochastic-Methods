[prob, stderr_of_est] = ObtuseTriangleProb(20000);

function [prob, stderr_of_est] = ObtuseTriangleProb(M)

%%% M samples of 3 points in the (x,y) plane %%%
    p = randn(M,2,3); 
    
%%% Calculate the three side lengths for all samples using
%%% the vector norm and orders each set of side lengths 
    lines = sort([vecnorm(p(:,:,1)-p(:,:,2),2,2), ...
        vecnorm(p(:,:,2)-p(:,:,3),2,2), ...
        vecnorm(p(:,:,3)-p(:,:,1),2,2)],2);

%%% Count the sets of points that meet the Pythagorean Thm. 
%%% of a^2 + b^2 < c^2 for obtuse triangles
    num_tri = sum(lines (:,1).^2 + lines(:,2).^2 < lines(:,3).^2);%0;

%%% Estimates probability
    prob = num_tri/M;
%%% Variance estimate of a Bernoulli Random Variable
    var_est = prob*(1-prob); 
%%% Standard error of estimate
    stderr_of_est = sqrt(var_est)/sqrt(M); 
     
end
function B_rand = signflip(B)
% Sign flip the values of B.

B_rand = B .* ((randi(2, size(B,1), size(B,2))-1).*2-1);
end
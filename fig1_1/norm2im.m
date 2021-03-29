function out = norm2im(vec,ref,n)

%Convert from normal to image coordinates

vec(vec < 0) = NaN; 

log_vec = log(vec); 
log_ref = log(ref);

adj1 = log_vec - min(log_ref);
adj2 = adj1/(max(log_ref) - min(log_ref));

out = adj2*(n-1) + 1;

end


% Low-Fidelity Model
function yL = low_fidelity(xL);
yL = 500.*xL.*(1-((xL.^2)/2));
end
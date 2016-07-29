% Calculate free energy based on Bins under certain temperature T. 
% Bins is the result from getBin

function F_energy = calFreeEnergy(Bins, T)
    Bins0 = Bins + 0.000001;
    F_energy = -T*log(Bins0)*0.001987; 
    F_energy = min(F_energy, 1);
end

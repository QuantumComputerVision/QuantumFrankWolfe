function [PisCells, energies, qubits] = load_dwave_text_output(fileName, N, NC)

D = dlmread(fileName);
cols = size(D,2);
bits = N^2*(NC-1);
I = eye(N);
I = I(:)';

energies = [];
qubits = [];
PisCells = {};
for i=1:2:size(D,1)
   remain = bits - cols;
   q = [I D(i,:) D(i+1,1:remain)]; 
   Pis = perms_q_to_cell(q, N);
   energy = D(i+1, remain+1);
   qubit =  D(i+1, remain+2);
   
   PisCells = [PisCells Pis];
   qubits = [qubits qubit];
   energies = [energies energy];
end

% sort by energy states
[energies, indSort] = sort(energies);

PisCells = PisCells(:, indSort);
qubits = qubits(indSort);
end

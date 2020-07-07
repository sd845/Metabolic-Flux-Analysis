using Pkg
Pkg.add("Convex")
Pkg.add("ECOS")
Pkg.add("SCS")
Pkg.add("Gurobi")

# Let us first make the Convex.jl module available
using Convex, ECOS, SCS, Gurobi
using DelimitedFiles

# reading stoichiometric matrix
S = readdlm("Network_deGFP.dat");
(S_rows,S_cols)=size(S)

#defining array of fluxes
v = [i for i in 1:S_cols]

#defining indices of measured fluxes: vm
vm = [j for j in 192:S_cols]

#defining indices of fluxes to be calculated: vc
vc=setdiff(v,vm);

#total number of flues
n = S_cols;
n_vc = length(vc);
n_vm = n - n_vc;

#size of S matrix
(S_rows,S_cols)=size(S)

#to get Sc matrix from S
Sc = zeros(S_rows,n_vc)
vm = sort(vm, rev=true)
for i in 1:S_rows
    A = S[i,:]
    for j in vm
        A = deleteat!(A,j)
    end
    Sc[i,:]=A
end
display(Sc);

n = S_cols;
n_vc = length(vc);
display(n_vc)
n_vm = n - n_vc;
#to get Sm matrix from S
Sm = zeros(S_rows,n_vm)
vc=setdiff(v,vm);
vc = sort(vc, rev=true)
display(vc)
for i in 1:S_rows
    A = S[i,:]
        for j in vc
        A = deleteat!(A,j)
    end
    Sm[i,:]=A
end
display(Sm);

#vm measured from experimental data
vm=readdlm("vm.dat");
display(vm)

#calculating b matrix
b=-Sm*vm;
display(b)
vc = Variable(n_vc);

# The problem is to minimize- subject to contraints
# This can be done by: minimize(objective, constraints)

#==============================================TXTL=====================================================#
RNAP_concentration_nM = 75;
RNAP_elongation_rate = 25;
RIBOSOME_concentration = 0.0016;
RIBOSOME_elongation_rate = 2;
kd = 5.2;
mRNA_length = 683;
protein_length = 229;
gene_copies = 3.125e10;
volume = 10e-6;
polysome_amplification = 10;
plasmid_saturation_coefficient = 3.5;
mRNA_saturation_coefficient = 0.045;
inducer = 35;
#====================================Transcription===================================================#
#Compute the promoter strength P -
hill_parameter = 1;
KD = 130;
K1 = 0.014;
K2 = 10;
K1_T7 = 10;
f = inducer^hill_parameter/(KD^hill_parameter+inducer^hill_parameter)
P = (K1+K2*f)/(1+K1+K2*f);
gene_concentration = gene_copies*(1e9/6.02e23)*(1/volume);
saturation_term = (gene_concentration)/(plasmid_saturation_coefficient+gene_concentration);
RNAP_concentration = RNAP_concentration_nM/1e6; #nM to mM
TX = (RNAP_elongation_rate*(1/mRNA_length)*(RNAP_concentration)*(saturation_term)*3600)*P;
#====================================Translation===================================================#
mRNA_steady_state = (TX/kd);
translation_rate_constant = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length)*3600;
TL = translation_rate_constant*RIBOSOME_concentration*mRNA_steady_state/(mRNA_saturation_coefficient+mRNA_steady_state);
#===================================================================================================#

problem = minimize(sumsquares(Sc*vc - b) + norm2(vc))
for i in 1:n_vc
    problem.constraints += vc[i] <= 100
    problem.constraints += vc[i] >= 0
end
problem.constraints += vc[167] == TX
problem.constraints += vc[169] == TX
problem.constraints += vc[170] <= TL
problem.constraints += Sm*vm + Sc*vc == 0

# Solve the problem by calling solve!
solve!(problem, SCSSolver(max_iter=30000))
#solve!(problem, ECOSSolver())

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval
println(round(problem.optval, digits=2))
println(round.(vc.value, digits=2))

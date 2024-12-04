import math
import operator
from qiskit import QuantumCircuit,QuantumRegister,ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit.compiler import transpile 
from qiskit.circuit.library import MCXGate 
from qiskit.visualization import *

from functions import to_binary
from functions import MaxCut



class FragmentedQGA:

    def __init__(self, graph, shots):
        self.graph = graph
        self.shots = shots
        self.individual_size = graph.number_of_nodes() #The total number of qubits needed for representing an individual
        self.fitness_size = 4 #math.ceil(math.log2(graph.number_of_edges())) #The number of qubits needed for representing the fitness value
        self.population_size = 2**self.individual_size #Superposition of all possible solutions
        self.grover_iterations = int(math.sqrt(2**(self.fitness_size))) #The maximum number of Grover iterations in which a solution should be found
        self.total_edges = graph.number_of_edges()

    def weight_check(self):
        """
        Checks if the graph has weighted edges and take them into account
        """
        has_unweighted_edges = any("weight" not in data for _, _, data in self.graph.edges(data=True))
        if has_unweighted_edges: 
            for u, v, data in self.graph.edges(data=True):
                if "weight" not in data:
                    data["weight"] = 1 # If unweighted, set weight as 1 by default

    def get_ufit_instruction(self):
        #define and initialize the individual quantum register
        ind_qreg = QuantumRegister(self.individual_size,"ind_qreg")
        #define and initialize the fitness quantum register. 
        fit_qreg = QuantumRegister(self.fitness_size+1,"fit_qreg")
        # Define the Multicontrolled Toffoli gates
        mcx_gate = MCXGate(num_ctrl_qubits=self.individual_size)
        #create the ufit subcircuit
        qc = QuantumCircuit(ind_qreg,fit_qreg,name="U$_fit$")
        for i in range(0,self.population_size):
            """
            For each individual in population, apply X-gate following the binary representation
            """
            individual_binary = to_binary(i, self.individual_size, True)
            for k in range(0,self.individual_size):
                if individual_binary[k] == 0:
                    qc.x(ind_qreg[k])
        
            partition = individual_binary
            #calculate cut value
            cut_value = MaxCut(self.graph).evaluate(partition)
            cut_value_bin = to_binary(cut_value,self.fitness_size,True)
    
            """
            Set the fitness value in fitness quantum register for each individual
            """
            for k in range(0,self.fitness_size):
                if cut_value_bin[k]==1:
                    qc.append(mcx_gate, [ind_qreg[j] for j in range(0, self.individual_size)] + [fit_qreg[k]]) #new
            # If the cut value is lower than the theoretical lower bound for the maxcut, ignore it.
            if cut_value >= self.total_edges // 2:
                qc.append(mcx_gate, [ind_qreg[j] for j in range(0, self.individual_size)] + [fit_qreg[self.fitness_size]]) #new
            #reset individual
            for k in range(0,self.individual_size):
                if individual_binary[k] == 0:
                    qc.x(ind_qreg[k])
            qc.barrier()
        return qc.to_instruction()

    def get_oracle_instruction(self, positive_value_array):
        #define and initialize fitness quantum register
        fit_reg = QuantumRegister(self.fitness_size,"fqreg")
        #define and initialize max quantum register
        no_of_edges_reg=QuantumRegister(self.fitness_size,"noqreg")
        #define and initialize carry quantum register
        carry_reg = QuantumRegister(3,"cqreg")
        #define and initialize oracle workspace quantum register
        oracle = QuantumRegister(1,"oqreg")
        #create Oracle subcircuit
        oracle_circ = QuantumCircuit(fit_reg,no_of_edges_reg,carry_reg,oracle,name="O")
        # Define the Multicontrolled Toffoli gates
        mcx_gate = MCXGate(num_ctrl_qubits=self.fitness_size)
        
        #define majority operator
        def majority(circ,a,b,c):
            circ.cx(c,b)
            circ.cx(c,a)
            circ.ccx(a, b, c)
        #define unmajority operator
        def unmaj(circ,a,b,c):
            circ.ccx(a, b, c)
            circ.cx(c, a)
            circ.cx(a, b)
        #define the Quantum Ripple Carry Adder
        def adder_4_qubits(p, a0, a1, a2, a3, b0, b1, b2, b3, cin, cout):
            majority(p, cin, b0, a0)
            majority(p, a0, b1, a1)
            majority(p, a1, b2, a2)
            majority(p, a2, b3, a3)
            p.cx(a3, cout)
            unmaj(p, a2, b3, a3)
            unmaj(p, a1, b2, a2)
            unmaj(p, a0, b1, a1)
            unmaj(p, cin, b0, a0)
        
        """
        Subtract max value. We start by storing the max value in the quantum register. Such, considering that 
        all qubits are |0>, if on position i in positive_value_array there's 0, then qubit i will be negated. Otherwise, 
        if on position i in positive_value_array there's a 1, by default will remain 0 in no_of_edges_reg quantum
        register. For performing subtraction, carry-in will be set to 1.
        """
        for i in range(0,self.fitness_size):
            if positive_value_array[i]==0:
                oracle_circ.x(no_of_edges_reg[i])
        oracle_circ.x(carry_reg[0])
    
        adder_4_qubits(oracle_circ, no_of_edges_reg[0],no_of_edges_reg[1],no_of_edges_reg[2],no_of_edges_reg[3],
                   fit_reg[0],fit_reg[1],fit_reg[2],fit_reg[3],
                   carry_reg[0],carry_reg[1]);
    
        
        oracle_circ.barrier()
        """
        Reset the value in no_of_edges_reg and carry-in
        """
        oracle_circ.x(no_of_edges_reg)
        oracle_circ.x(carry_reg[0])
        
        """
        Mark the corresponding basis states by shifting their amplitudes.
        """
        
        oracle_circ.h(oracle[0])
        oracle_circ.append(mcx_gate, [fit_reg[j] for j in range(0, self.fitness_size)] + [oracle[0]]) #new
        oracle_circ.h(oracle[0])
        
        """
        Restore the fitness value by adding max value.
        """
        adder_4_qubits(oracle_circ, no_of_edges_reg[0],no_of_edges_reg[1],no_of_edges_reg[2],no_of_edges_reg[3],
                   fit_reg[0],fit_reg[1],fit_reg[2],fit_reg[3],
                   carry_reg[0],carry_reg[2]);
        return oracle_circ.to_instruction()

    def get_grover_iteration_subcircuit(self):
        #define and initialize fitness quantum register
        fit_qreg = QuantumRegister(self.fitness_size+1,"fqreg")
        #define and initialize oracle workspace quantum register
        oracle_ws = QuantumRegister(1,"ows")
        #create grover diffuser subcircuit
        grover_circ = QuantumCircuit(fit_qreg,oracle_ws,name ="U$_s$")
        # Define the Multicontrolled Toffoli gates
        mcx_gate = MCXGate(num_ctrl_qubits=self.fitness_size + 1)
    
        grover_circ.h(fit_qreg)
        grover_circ.x(fit_qreg)
    
        grover_circ.h(oracle_ws[0])
    
        grover_circ.append(mcx_gate, list(range(self.fitness_size + 1)) + [oracle_ws[0]]) #new
    
        grover_circ.h(oracle_ws[0])
    
    
        grover_circ.x(fit_qreg)
        grover_circ.h(fit_qreg)
        grover_circ.h(oracle_ws)
    
        return grover_circ.to_instruction()


    def run(self):
        #calculate the number of edges in graph
        self.weight_check()
       # print("Number of vertices:{0}".format(self.individual_size))
       # print("Number of edges:{0}".format(self.total_edges))
       # print("Number of Grover Iterations:{0}".format(self.grover_iterations))
       # print("Population Size in Superposition:{0}".format(self.population_size))
        #define a list for storing the results
        final_results = []
        
        #Start with one Grover iteration and run the algorithm with different number of iterations
        for iterations in range(1,self.grover_iterations+1):
            # print("Iteration #{0}".format(iterations))
            
           # print("Preparing quantum registers and creating quantum circuit...")
            ind_qreg = QuantumRegister(self.individual_size,"ireg")
            fit_qreg = QuantumRegister(self.fitness_size+1,"freg") #8 qubits fitness + 1 valid
            carry_qreg = QuantumRegister(2*self.grover_iterations+1,"qcarry")
            oracle = QuantumRegister(1,"oracle")
            creg = ClassicalRegister(self.individual_size,"reg")
            no_of_edges_qreg = QuantumRegister(self.fitness_size,"pos_max_qreg")
            qc = QuantumCircuit(ind_qreg,fit_qreg,carry_qreg,oracle,no_of_edges_qreg,creg)

           # print("Setting Individuals in Superposition...")
            qc.h(ind_qreg)
            qc.h(oracle)
    
            pos_value_bin = to_binary(self.total_edges, self.fitness_size,True)
    
    
            ufit_instr = self.get_ufit_instruction()
            oracle_instr = self.get_oracle_instruction(pos_value_bin)
            grover_iter_inst = self.get_grover_iteration_subcircuit()
    
            qc.append(ufit_instr, [ind_qreg[q] for q in range(0,self.individual_size)]+
                                  [fit_qreg[q] for q in range(0,self.fitness_size+1)]
                    )
            
            for it in range(0,iterations):
                
                qc.append(oracle_instr,[fit_qreg[q] for q in range(0,self.fitness_size)]+
                                   [no_of_edges_qreg[q] for q in range(0,self.fitness_size)]+
                                   [carry_qreg[0],carry_qreg[2*it+1],carry_qreg[2*it+2]]+
                                   [oracle[0]])
                qc.append(grover_iter_inst, [fit_qreg[q] for q in range(0,self.fitness_size+1)]+[oracle[0]])
            
           # print("Setting up Circuits...")
            qc.measure(ind_qreg,creg)
            
            simulation_results = []
            
            backend = AerSimulator()
    
           # print("Starting simulator...")
            # Perform 10 measurements for each circuit
            for run in range(0,10):
                try:
                    mapped_circuit = transpile(qc, backend=backend)
                    job = backend.run(mapped_circuit, shots=self.shots)
                    job.status()
                    #job_monitor(job)
                    results = job.result()
                    answer = results.get_counts()
                    #Get the result with the maximum number of counts
                    max_item =max(answer.items(), key=operator.itemgetter(1))
                   # print(max_item)
                    #Store the result.
                    simulation_results.append(max_item)
                except Exception as e:
                    # print(str(e))
                    print("Error on run {0} with {1} grover iterations".format(run,iterations))
            final_results.append((iterations,simulation_results))
        return final_results


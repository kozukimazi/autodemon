import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.linalg import logm 
from scipy import integrate
import cmath
import os
from scipy.sparse import kron, eye, csr_matrix


class FermionicSystem:
    def __init__(self, n_sites=3):
        self.n_sites = n_sites
        self.dim = 2**n_sites
        self.sigma_z = csr_matrix([[1, 0], [0, -1]])
        self.sigma_plus = csr_matrix([[0, 1], [0, 0]])  # Creation
        self.sigma_minus = csr_matrix([[0, 0], [1, 0]]) # Annihilation
        
        # Precompute Jordan-Wigner operators
        self.creation_ops = [self.build_operator(i, is_creation=True) for i in range(n_sites)]
        self.annihilation_ops = [self.build_operator(i, is_creation=False) for i in range(n_sites)]
        
        # Build Fock basis states
        self.fock_basis = self.build_fock_basis()
    
    def build_operator(self, site, is_creation=True):
        """
        Build fermionic operator for given site using Jordan-Wigner transformation
        """
        op = csr_matrix((1,1))  # Initialize with 1x1 identity
        
        # Build the string of sigma_z operators for sites before the target site
        for i in range(site):
            op = kron(op, self.sigma_z, format='csr')
        
        # Add the creation/annihilation operator at the target site
        if is_creation:
            op = kron(op, self.sigma_plus, format='csr')
        else:
            op = kron(op, self.sigma_minus, format='csr')
        
        # Fill remaining sites with identity
        for i in range(site+1, self.n_sites):
            op = kron(op, eye(2), format='csr')
        
        return op
    
    def build_fock_basis(self):
        """Build all possible occupation number states"""
        basis = []
        for i in range(self.dim):
            # Convert index to binary representation
            state = [int(x) for x in format(i, f'0{self.n_sites}b')]
            basis.append(state)
        return basis
    
    def get_state_vector(self, occupation):
        """
        Convert occupation list to state vector
        e.g., [1,0,1] → vector for |1⟩⊗|0⟩⊗|1⟩
        """
        if len(occupation) != self.n_sites:
            raise ValueError("Occupation list length must match number of sites")
            
        index = int(''.join(map(str, occupation)), 2)
        psi = np.zeros(self.dim)
        psi[index] = 1
        return psi
    
    def partial_trace(self, rho, keep_sites):
        """
        Partial trace with proper fermionic sign handling
        rho: density matrix (2^n × 2^n)
        keep_sites: list of sites to keep (e.g., [0,1])
        """
        trace_sites = [i for i in range(self.n_sites) if i not in keep_sites]
        n_keep = len(keep_sites)
        dim_keep = 2**n_keep
        
        # Initialize reduced density matrix
        rho_reduced = np.zeros((dim_keep, dim_keep), dtype=complex)
        
        # Iterate over all basis states of the kept subsystem
        for i_keep in range(dim_keep):
            for j_keep in range(dim_keep):
                total = 0
                
                # Iterate over all possible states of traced-out sites
                for trace_state in range(2**len(trace_sites)):
                    # Construct full basis indices
                    full_i = self.reconstruct_full_index(i_keep, trace_state, keep_sites, trace_sites)
                    full_j = self.reconstruct_full_index(j_keep, trace_state, keep_sites, trace_sites)
                    
                    # Get sign factor from fermionic reordering
                    sign = self.get_reordering_sign(full_i, full_j, keep_sites, trace_sites)
                    
                    total += rho[full_j, full_i] * sign
                
                rho_reduced[j_keep, i_keep] = total
        
        return rho_reduced
    
    def reconstruct_full_index(self, keep_index, trace_index, keep_sites, trace_sites):
        """Reconstruct the full basis index from subsystem indices"""
        full_state = [0]*self.n_sites
        
        # Fill in kept sites
        keep_bits = format(keep_index, f'0{len(keep_sites)}b')
        for pos, bit in zip(keep_sites, keep_bits):
            full_state[pos] = int(bit)
            
        # Fill in traced-out sites
        trace_bits = format(trace_index, f'0{len(trace_sites)}b')
        for pos, bit in zip(trace_sites, trace_bits):
            full_state[pos] = int(bit)
            
        return int(''.join(map(str, full_state)), 2)
    
    def get_reordering_sign(self, i, j, keep_sites, trace_sites):
        """
        Compute the sign factor from fermionic reordering when performing partial trace
        """
        # Get occupation numbers for both states
        occ_i = [int(x) for x in format(i, f'0{self.n_sites}b')]
        occ_j = [int(x) for x in format(j, f'0{self.n_sites}b')]
        
        sign = 1
        
        # For each traced-out site, count how many kept fermions are to its left
        for t in trace_sites:
            for k in keep_sites:
                if k < t and (occ_i[k] == 1 or occ_j[k] == 1):
                    sign *= -1
                    
        return sign
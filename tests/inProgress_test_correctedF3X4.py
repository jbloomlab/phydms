import scipy
import scipy.optimize
import unittest
from phydmslib.constants import *

class testYNGKP_M0_correctedF3X4(unittest.TestCase):

    def testYNGKP_M0_correctedF3X4(self):
        """Test `scipy.optimize.root` function fo the corrected F3X4 frequencies."""

        scipy.random.seed(1)
        self.e_pw = scipy.random.uniform(0.12, 1.0, size = (3, N_NT))
        self.e_pw = self.e_pw/self.e_pw.sum(axis=1, keepdims=True)

        def F(phi_pw):
            phi_pw_reshape = phi_pw.reshape((3,4)) #reshape flattened array
            functionList = [] #the list of 12 equations. returned at the end
            stop_frequency = [] #frequency of each of the three stop codons based on the phi_pw values

            for x in range(N_STOP):
                codonFrequency = STOP_CODON_TO_NT_INDICES[x] * phi_pw_reshape
                codonFrequency = scipy.prod(codonFrequency.sum(axis=1)) #sum across each position (only one value per row will be non-zero)
                stop_frequency.append(codonFrequency) #multiply the phi_pw frequency for each nt in codon
            C = scipy.sum(stop_frequency) # this is the denominator from the Pond paper

            for p in range(3):
                for w in range(N_NT):
                    s = 0
                    for x in range(N_STOP):
                        if STOP_CODON_TO_NT_INDICES[x][p][w] == 1: #if the nucleotide, positition pair is in the codon
                            s += stop_frequency[x] #add the frequency of the codon
                    functionList.append((phi_pw_reshape[p][w] - s)/(1 - C) - self.e_pw[p][w]) #final calc
            return functionList

        phi_pw = self.e_pw.copy().flatten()
        assert((phi_pw.reshape((3,N_NT)) == self.e_pw).all())
        with scipy.errstate(invalid='ignore'):
            result = scipy.optimize.root(F, phi_pw,
                    tol=1e-8)
            assert result.success, "Failed: {0}".format(result)
            phi_pw = result.x.reshape((3,N_NT))
            self.assertTrue(scipy.allclose(phi_pw.sum(axis = 1),
                scipy.ones(3, dtype='float'),atol=1e-4, rtol=5e-3))
            assert(phi_pw[0][3] >= self.e_pw[0][3]) #The calculated frequency for T at position 0 should be higher than the empirical frequency
            assert(phi_pw[0][0] <= self.e_pw[0][0]) #The calculated frequency for A at position 0 should be lower than the numerical frequency

            newValues = F(phi_pw.flatten())
            self.assertTrue(scipy.allclose(scipy.asarray(newValues).reshape((3,4)),
                scipy.zeros((3, N_NT), dtype='float'), atol=1e-4, rtol=5e-3))




if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)

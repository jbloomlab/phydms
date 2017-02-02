import scipy
import scipy.optimize
import unittest
from phydmslib.constants import *

class testYNGKP_M0_correctedF3X4(unittest.TestCase):

    def testYNGKP_M0_correctedF3X4(self):
        """Initialize `YNGKP_M0`, test values, update, test again."""

        scipy.random.seed(1)
        self.e_pw = scipy.random.uniform(0.12, 1.0, size = (3, N_NT))
        self.e_pw = self.e_pw/self.e_pw.sum(axis=1, keepdims=True)

        def F(phi_pw):
            #the index of phi_x for pos `p` and nt `w` is (4*p + w) becaue of the flattend array
            functionList = []
            pi_x = 3 * phi_pw[0 * 4 + 3] + 2 * phi_pw[1 * 4 + 0] + phi_pw[1 * 4 + 2] + 2 * phi_pw[2 * 4 + 0] + phi_pw[2 * 4 + 2] #part of the denom
            for p in range(3): #for each position
                for w in range(4): #for each nucleotide, #this matches the flatten/reshape order
                    s = 0
                    for stop in range(N_STOP): #for each stop codon
                        if(NT_CODON_POS_TO_STOP[p][w][stop]): #if `pw` part of stop codon
                            ntList = [NT_TO_INDEX[STOP[stop][y]] for y in range(3)] #make a list (in codon position order) of the nt indexes of that stop codon
                            freqList = scipy.asarray([phi_pw[y * 4 + ntList[y]] for y in range(3)]) #get the phi values
                            s += (scipy.prod(freqList)) #multiply the phi values and add to s
                    functionList.append(((phi_pw[p * 4 + w] - s / (1 - pi_x)) - self.e_pw[p][w])) #do the final calculations
            return(scipy.asarray(functionList))

        print(self.e_pw)
        phi_pw = self.e_pw.copy().flatten()
        assert((phi_pw.reshape((3,4)) == self.e_pw).all())
        with scipy.errstate(invalid='ignore'):
            result = scipy.optimize.root(F, phi_pw,
                    tol=1e-8)
            assert result.success, "Failed: {0}".format(result)
            phi_pw = result.x.reshape((3,4))
        print()
        print(phi_pw)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)

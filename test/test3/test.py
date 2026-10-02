import numpy as np 

#reference_data = np.loadtxt('plt_level_convergence_6000_3e4_copy.dat')
reference_data = np.loadtxt('plt_level_convergence_test.dat')

test_data = np.loadtxt('plt_level_convergence.dat')

tolerance = 3e-2

def testArray(baselineArray, testdataArray, testname):

    #select = baselineArray > 1e-2 * np.max(baselineArray) 
    
    percentDiff = np.abs( ((baselineArray - testdataArray) / baselineArray) )

    if np.any(percentDiff > tolerance):
        print(f'{testname} test failed.',np.max(percentDiff))
        for jj in range(0,len(percentDiff)):
            print('     Entry {:2}. Calc: {:10.6e} Ref: {:10.6e}'.format(jj,testdataArray[jj],baselineArray[jj]))
            

    else:
        print(f'{testname} test passed.',np.max(percentDiff))


testArray(reference_data[:,1],test_data[:,1],'PLT Convergence')


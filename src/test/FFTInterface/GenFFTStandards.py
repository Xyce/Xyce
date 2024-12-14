
#
# This is a python program using Numpy's FFT routines
# to calculate inputs and outputs to test Xyce's FFTInterface
#
# Signals are generated and saved to CSV file for FFTUnitTests to use
#
# The signals are:
# 1. Single frequency and multi frequency periodic functions.
# 2. Non periodic functions such as gaussian and step functions.
# 3. Have a power of 2 number of points
# 4. Have a Non-power of 2 number of points 
#
# Data is generated in the time domain, transformed to the frequency domain,
# and back transformed to the time domain.  Data from all three stages is 
# saved to a CSV file

import math
import numpy as np
import matplotlib.pyplot as pt

def SingleFrequencySin( timeVec, periods=3.25, amplitude=1.5, phase=(math.pi/10) ):
  outputVec = np.zeros_like(timeVec)
  deltaTime = max(timeVec)-min(timeVec)
  frequency = 2.0*math.pi/deltaTime;
  outputVec = amplitude*np.sin( periods*frequency*timeVec + phase)
  return outputVec
    
def TwoFrequencySin( timeVec, periods=3.25, amplitude=1.5, phase=(math.pi/10) ):
  outputVec = np.zeros_like(timeVec)
  deltaTime = max(timeVec)-min(timeVec)
  frequency = 2.0*math.pi/deltaTime;
  outputVec = amplitude*np.sin( periods*frequency*timeVec + phase) + 0.5*amplitude*np.sin( periods*4*frequency*timeVec + phase)
  return outputVec
  
def ThreeFrequencySin( timeVec, periods=3.25, amplitude=1.5, phase=(math.pi/10) ):
  outputVec = np.zeros_like(timeVec)
  deltaTime = max(timeVec)-min(timeVec)
  frequency = 2.0*math.pi/deltaTime;
  outputVec = amplitude*np.sin( periods*frequency*timeVec + phase) + 0.5*amplitude*np.sin( periods*4*frequency*timeVec + 1.1*phase) + 0.25*amplitude*np.sin( periods*7*frequency*timeVec + 2.5*phase)
  return outputVec
  
def GaussianShape( timeVec,  E01=10.0,  a1=2e4, b1=3.0e5 ):
  outputVec = np.zeros_like(timeVec)
  outputVec = E01 * ( np.exp( -a1*timeVec ) - np.exp( -b1*timeVec ) )
  return outputVec
  
def SingleFreqStep( timeVec, periods=3.25, amplitude=1.5, phase=(math.pi/10) ):
  outputVec = np.zeros_like(timeVec)
  deltaTime = max(timeVec)-min(timeVec)
  frequency = 2.0*math.pi/deltaTime;
  outputVec = np.clip(amplitude*np.sin( periods*frequency*timeVec + phase), 0.0, 0.5)
  return outputVec
  
def DoCalculationsAndSave(FileName, timeVec, functionVec, showplots=False ):
  functionFFT = np.fft.fft(functionVec)
  functionIFFT = np.fft.ifft(functionFFT)
  np.savetxt(FileName, (timeVec, functionVec, np.real(functionFFT), np.imag(functionFFT), np.real(functionIFFT), np.imag(functionIFFT)))
  if( showplots ):
    pt.plot(timeVec, functionVec)
    pt.title('Original Function')
    pt.show()
    pt.plot(np.real(functionFFT))
    pt.plot(np.imag(functionFFT))
    pt.title('FFT Function')
    pt.show()
    pt.plot(timeVec, functionVec)
    pt.plot(timeVec, np.real(functionIFFT), 'o')
    pt.title('Original and IFFT Function')
    pt.show()
    
def CalculateFFTStandards():
  print("In CalculateFFTStandards")
  numPointsPowerOfTwo = pow(2,8)
  numPointsNonPowerOfTwo = numPointsPowerOfTwo + 17
  timeBasePowerOfTwo = np.linspace(0, 5e-3, num=numPointsPowerOfTwo)
  timeBaseNonPowerOfTwo = np.linspace(0, 5e-3, num=numPointsNonPowerOfTwo)

  periodicOneFreqPowerOfTwo = SingleFrequencySin(timeBasePowerOfTwo, amplitude=5.0, phase=(math.pi/33))
  DoCalculationsAndSave('Sin1fEven.txt', timeBasePowerOfTwo, periodicOneFreqPowerOfTwo)
  
  periodicOneFreqNonPowerOfTwo = SingleFrequencySin(timeBaseNonPowerOfTwo, amplitude=4.0, phase=(math.pi/3))
  DoCalculationsAndSave('Sin1fOdd.txt', timeBaseNonPowerOfTwo, periodicOneFreqNonPowerOfTwo)
    
  periodicTwoFreqPowerOfTwo = TwoFrequencySin(timeBasePowerOfTwo, amplitude=5.0, phase=(math.pi/33))
  DoCalculationsAndSave('Sin2fEven.txt', timeBasePowerOfTwo, periodicTwoFreqPowerOfTwo)
  
  periodicTwoFreqNonPowerOfTwo = TwoFrequencySin(timeBaseNonPowerOfTwo, amplitude=4.0, phase=(math.pi/3))
  DoCalculationsAndSave('Sin2fOdd.txt', timeBaseNonPowerOfTwo, periodicOneFreqNonPowerOfTwo)
  
  periodicThreeFreqPowerOfTwo = ThreeFrequencySin(timeBasePowerOfTwo, amplitude=5.0, phase=(math.pi/33))
  DoCalculationsAndSave('Sin3fEven.txt', timeBasePowerOfTwo, periodicThreeFreqPowerOfTwo)
  
  periodicThreeFreqNonPowerOfTwo = TwoFrequencySin(timeBaseNonPowerOfTwo, amplitude=4.0, phase=(math.pi/3))
  DoCalculationsAndSave('Sin3fOdd.txt', timeBaseNonPowerOfTwo, periodicThreeFreqNonPowerOfTwo)
  
  gaussianPowerOfTwo = GaussianShape(timeBasePowerOfTwo)
  DoCalculationsAndSave('Gauss1Even.txt', timeBasePowerOfTwo, gaussianPowerOfTwo)
  
  gaussianNonPowerOfTwo = GaussianShape(timeBaseNonPowerOfTwo)
  DoCalculationsAndSave('Gauss1Odd.txt', timeBaseNonPowerOfTwo, gaussianNonPowerOfTwo)
  
  gaussianLongPowerOfTwo = GaussianShape(timeBasePowerOfTwo,a1=2e3, b1=3.0e4)
  DoCalculationsAndSave('Gauss2Even.txt', timeBasePowerOfTwo, gaussianLongPowerOfTwo, showplots=False)
  
  gaussianLongNonPowerOfTwo = GaussianShape(timeBaseNonPowerOfTwo,a1=2e3, b1=3.0e4)
  DoCalculationsAndSave('Gauss2Odd.txt', timeBaseNonPowerOfTwo, gaussianLongNonPowerOfTwo)

  squarePowerOfTwo = SingleFreqStep(timeBasePowerOfTwo)
  DoCalculationsAndSave('Step1Even.txt', timeBasePowerOfTwo, squarePowerOfTwo)
  squareNonPowerOfTwo = SingleFreqStep(timeBaseNonPowerOfTwo)
  DoCalculationsAndSave('Step1Odd.txt', timeBaseNonPowerOfTwo, squareNonPowerOfTwo)
 
#
# if this file is called directly then run the tests
#
if __name__ == '__main__':
  retval = CalculateFFTStandards()
  exit( retval)
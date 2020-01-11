"""
data in filename type:  v  I
Ignore the first row
by  mhf 20191103

"""


import numpy as np
import math



filename='inputdata.txt'

# def get_morse_frequency(frequency,vibration_intensity,excited_state_num,temperture):
#     # t  = temperture
#     # I  = vibration_intensity
#     n  = excited_state_num
#     h  = 6.626e(-34)
#     kb = 1.38e(-23)
#     v  = frequency*1e12   #real frequecy
#     k  = -h*v/(kb*temperture)
#     f  = exp(k*(0.5+0.98*n))
#

temperature = 80
excited_state_num = 10
anharmonicity_coefficient = 0.04
t  = temperature
n  = excited_state_num
p  = anharmonicity_coefficient
h  = 6.626e-34
kb = 1.38e-23

"""
Generate a series of frequecies with certain coefficients.
Frequency's dimensionality at certain frequency equals  excited_state_num

"""

def new_v(frequency,anharmonicity_coefficient,excited_state_num):
    p=anharmonicity_coefficient
    frequency=np.array(frequency)
    v1=[]
    for i in range(excited_state_num):
        v1.append((1-i*p)*frequency)
    return v1


"""
Generate a series of intensities with certain coefficients.
Intensities's dimensionality at certain frequency equals  excited_state_num

"""
def new_I(frequency,intensity,temperature,anharmonicity_coefficient,excited_state_num):
    p = anharmonicity_coefficient
    k = -h * frequency / (kb * temperature)  # coefficient of f_bolzman function's exponent
    intensity=np.array(intensity)
    new_f_b=[]     #bolzman function of this
    for i in range(excited_state_num):
        f_b=np.exp(k*(0.5+(1-p)*i))
        new_f_b.append(f_b)
    pop_prop=[]
    sum_f_b = np.sum(new_f_b)
    for i in range(excited_state_num):
         pop_prop.append(new_f_b[i]/sum_f_b)
    I1=intensity*pop_prop
    return I1




vs,Is=np.loadtxt(filename, skiprows=1,unpack=True)  #vs:frequencies, Is: intencities
vs  = vs*1e12   #real frequecy
#k  = -h*vs/(kb*t)   #coefficient of f_bolzman function
#f_bs = np.exp(k*(0.5+(1-p)*n))   #function of bolzman


v1 = []; I1=[]
for (v, I) in zip(vs, Is):

    v1.append(new_v(v, p, n))
    I1.append(new_I(v,I,t,p,n))

v1=np.array(v1).reshape(np.size(v1),1)   #new morse frequency
I1=np.array(I1).reshape(np.size(I1),1)   #new morse intensity

v_I_array=np.concatenate((v1/1e12, I1),axis=1)


np.savetxt('result.txt', v_I_array, fmt='%.3f')

# pop_prop=[]
# sum_f_b=np.sum(f_b)
# for i in range(3):
#     pop_prop.append(f_b[i]/sum_f_b)        #poulation propertion




import numpy as np
import matplotlib.pyplot as plt

def get_correlation_from_value(ref_data, target_data, value=0.0):
    correlated_data_sup = []
    correlated_data_inf = []

    for ref, data in zip(ref_data, target_data):
        if ref < value:
            correlated_data_sup.append(data)
        else:
            correlated_data_inf.append(data)

    return correlated_data_inf, correlated_data_sup


slides = np.loadtxt('slides.dat')

# LONGITUDINAL
inf, sup = get_correlation_from_value(slides[:, 1], slides[:, 0], value=2)
#np.savetxt('slides_inf_x.dat', inf)
#np.savetxt('slides_sup_x.dat', sup)
#plt.plot(inf)
#plt.plot(sup)
#plt.show()
plt.title('longitudinal')
plt.hist(slides[:, 0], histtype=u'step', label='full')
plt.hist(inf,histtype=u'step', label='lateral < 2')
plt.hist(sup,histtype=u'step', label='lateral > 2')
plt.legend()
plt.show()


# LATERAL
inf, sup = get_correlation_from_value(slides[:, 0], slides[:, 1], value=4)
#np.savetxt('slides_inf_y.dat', inf)
#np.savetxt('slides_sup_y.dat', sup)
#plt.plot(inf)
#plt.plot(sup)
#plt.show()
plt.title('lateral')
plt.hist(slides[:, 1], histtype=u'step', label='full')
plt.hist(inf, histtype=u'step', label='longitudinal < 4')
plt.hist(sup,histtype=u'step', label='longitudinal > 4')
plt.legend()
plt.show()


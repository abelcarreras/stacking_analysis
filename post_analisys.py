import numpy as np
import matplotlib.pyplot as plt


def get_correlation_from_value(ref_data, target_data, value=0.0):
    correlated_data_sup = []
    correlated_data_inf = []

    for ref, data in zip(ref_data, target_data):
        if ref < value:
            correlated_data_inf.append(data)
        else:
            correlated_data_sup.append(data)

    return correlated_data_inf, correlated_data_sup


dir = ''
slides = np.loadtxt(dir + 'slides.dat')
rotation_angle = np.loadtxt(dir + 'rotation.dat')
distance = np.loadtxt(dir + 'distance.dat')
neighbors = np.loadtxt(dir + 'neighbors.dat')


# LONGITUDINAL
inf, sup = get_correlation_from_value(slides[:, 1], slides[:, 0], value=6.5)
#np.savetxt('slides_inf_x.dat', inf)
#np.savetxt('slides_sup_x.dat', sup)
#plt.plot(inf)
#plt.plot(sup)
#plt.show()
plt.title('longitudinal')
N, bins, patches = plt.hist(slides[:, 0], histtype=u'step', label='full', normed=1)
plt.hist(inf, bins=bins, histtype=u'step', label='lateral < 5', normed=1)
plt.hist(sup, bins=bins, histtype=u'step', label='lateral > 5', normed=1)
plt.legend()
plt.show()


# LATERAL
inf, sup = get_correlation_from_value(slides[:, 1], slides[:, 1], value=6.5)
#np.savetxt('slides_inf_y.dat', inf)
#np.savetxt('slides_sup_y.dat', sup)
#plt.plot(inf)
#plt.plot(sup)
#plt.show()
plt.title('lateral')
N, bins, patches = plt.hist(slides[:, 1], histtype=u'step', label='full', normed=1)
plt.hist(inf, bins=bins, histtype=u'step', label='longitudinal < 4',normed=1)
plt.hist(sup, bins=bins, histtype=u'step', label='longitudinal > 4',normed=1)
plt.legend()
plt.show()



# Plot data

plt.title('Rotation angle')
plt.ylim([0,90])
plt.plot(rotation_angle, 'o')
plt.show()

inf, sup = get_correlation_from_value(slides[:, 1], rotation_angle, value=6.5)
plt.title('Rotation angle')
N, bins, patches = plt.hist(rotation_angle, histtype=u'step', bins=10, label='full', normed=1)
plt.hist(inf, bins=bins, histtype=u'step', label='lateral < 2', normed=1)
plt.hist(sup, bins=bins, histtype=u'step', label='lateral > 2', normed=1)
plt.legend()
plt.show()


plt.title('Distance')
plt.ylim([0, 10])
plt.plot(distance, 'o')
plt.show()

#inf, sup = get_correlation_from_value(np.linalg.norm(slides, axis=1), distance, value=2)
inf, sup = get_correlation_from_value(slides[:, 1], distance, value=6.5)

plt.title('Distance')
N, bins, patches = plt.hist(distance, histtype=u'step', bins=20, normed=1, label='full')
plt.hist(inf, bins=bins, histtype=u'step', normed=1,  label='displacement < 2')
plt.hist(sup, bins=bins, histtype=u'step', normed=1,  label='displacement > 2')
plt.legend()
plt.show()


plt.title('neighbors')
plt.plot(neighbors)
plt.show()

plt.title('neighbors')
plt.hist(neighbors, normed=1)
#plt.xticks(np.arange(-0.5, 8.5, 1.0), [str(k-1) for k in range(20)])
plt.show()

#plt.hist(rotation_angle[500:], normed=1)


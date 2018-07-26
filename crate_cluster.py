import numpy as np

# create cluster

# read structure


def read_txyz(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()

    symbols = []
    coordinates = []
    types = []
    connectivity = []

    for line in lines[1:]:
        data_list = line.split()

        symbols.append(data_list[1])
        coordinates.append(data_list[2:5])
        types.append(data_list[5])
        connectivity.append([int(val) for val in data_list[6:]])
    coordinates = np.array(coordinates, dtype=float)
    types = np.array(types, dtype=int)

    return symbols, coordinates, types, connectivity


def write_txyz(filename, symbols, coordinates, types, connectivity):

    f = open(filename, 'w')
    n_atoms = len(symbols)

    f.write('{}\n'.format(n_atoms))
    for i in range(n_atoms):
        f.write(' {} {} '.format(i+1, symbols[i]) +
                ' {0:12.6f} {1:12.6f} {2:12.6f} '.format(*coordinates[i]) +
                ' {} '.format(types[i]) +
                ' '.join([ str(val) for val in connectivity[i]]) +
                '\n')

    return symbols, coordinates, types, connectivity


def rotation_matrix(phi):
    return np.array([[np.cos(phi), -np.sin(phi)],
                     [np.sin(phi), np.cos(phi)]])

def mod_coordinates(symbols, coordinates, types, connectivity):

    n_atoms = len(symbols)
    s = symbols
    t = types
    ic = coordinates
    ci = list(connectivity)

    index = 0
    for i3 in range(1): #4
        for i2 in range(1):  #4
            for i in range(2):  #8
                if i3 == 0 and i2 == 0 and i == 0:
                    continue
                c = ic + np.array([25.0*(i3+0)+ i*4 , 25.0*(i2+0)+ i*5, 10.0*(i+0)])
                #print(np.array([0.0, 0.0, 10.0*(i+1.0)]))
                angle = np.random.random_sample()*np.pi*2.0
                angle = (i+1) * np.pi * 2 * 13.0 / 360.0
                #angle = np.pi * 2 * 70.0 / 360.0
                angle = 0.0
                c[:, :2] = np.dot(c[:, :2], rotation_matrix(angle))

                index += 1
                cn = []
                for j, k in enumerate(ci):
                    #cn.append((np.array(k) + n_atoms*(i+1)).tolist())
                    cn.append((np.array(k) + n_atoms*(index)).tolist())

                symbols = np.hstack([symbols, s])
                types = np.hstack([types, t])
                connectivity = np.hstack([connectivity, cn])
                coordinates = np.vstack([coordinates, c])

    connectivity = list(connectivity.tolist())

    return symbols, coordinates, types, connectivity


#symbols, coordinates, types, connectivity = read_txyz('qt4c_R.txyz')
symbols, coordinates, types, connectivity = read_txyz('QT4C/qt4c_R.txyz')

symbols, coordinates, types, connectivity = mod_coordinates(symbols, coordinates, types, connectivity)

#write_txyz('qt4c_helix.txyz', symbols, coordinates, types, connectivity)
write_txyz('QT4C/qt4c_test.txyz', symbols, coordinates, types, connectivity)


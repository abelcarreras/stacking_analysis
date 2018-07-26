import numpy as np


def reading_from_arc_file(file_name, nrelax=0, nmax=19000):

    tinker_file = open(file_name, 'r')

    trajectory = []

    number_of_atoms = int(tinker_file.readline().split()[0])
    # fpos = tinker_file.tell()
    check = tinker_file.readline().split()[0]
    if check == '1':
        is_cell = False
    else:
        is_cell = True

    tinker_file.seek(0)

    eof = False
    n=0
    while not eof:

        tinker_file.readline()
        if is_cell:
            tinker_file.readline()

        try:
            atomic_numbers = []
            atomic_elements = []
            atom_types = []
            connectivity = []
            coordinates = []
            for i in range(number_of_atoms):
                line = tinker_file.readline().split()
                coordinates.append(line[2:5])
                atomic_numbers.append(int(line[0]))
                atomic_elements.append(line[1])
                atom_types.append(line[5])
                connectivity.append([int(f) for f in line[6:]])

            if n >= nrelax:
                trajectory.append(np.array(coordinates, dtype=float))

        except IndexError:
            eof = True

        if n > nmax:
            break
        else:
            n += 1

    print(len(trajectory))
    tinker_file.close()
    return {'trajectory': np.array(trajectory),
            'atom_types': np.array(atom_types, dtype=int)[None].T,
            'atomic_elements': np.array(atomic_elements, dtype=str)[None].T,
            'atomic_numbers': np.array(atomic_numbers, dtype=int)[None].T,
            'connectivity': list(connectivity)}


def reading_from_arc_file_mmap(file_name, nrelax=0, nmax=19000):

    import mmap

    with open(file_name, 'r+') as tinker_file:

        file_map = mmap.mmap(tinker_file.fileno(), 0)

        trajectory = []

        number_of_atoms = int(file_map.readline().split()[0])
        #fpos = file_map.tell()
        check = file_map.readline().split()[0]
        if check == '1':
            is_cell = False
        else:
            is_cell = True

        file_map.seek(0)

        eof = False
        n = 0
        while not eof:

            file_map.readline()
            if is_cell:
                file_map.readline()

            try:
                atomic_numbers = []
                atomic_elements = []
                atom_types = []
                connectivity = []
                coordinates = []
                for i in range(number_of_atoms):
                    line = file_map.readline().split()
                    coordinates.append(line[2:5])
                    atomic_numbers.append(int(line[0]))
                    atomic_elements.append(line[1])
                    atom_types.append(line[5])
                    connectivity.append([int(f) for f in line[6:]])

                if n >= nrelax:
                    trajectory.append(np.array(coordinates, dtype=float))

            except IndexError:
                eof = True

            if n > nmax:
                break
            else:
                n += 1

        print(len(trajectory))
        file_map.close()

    return {'trajectory': np.array(trajectory),
            'atom_types': np.array(atom_types, dtype=int)[None].T,
            'atomic_elements': np.array(atomic_elements, dtype=str)[None].T,
            'atomic_numbers': np.array(atomic_numbers, dtype=int)[None].T,
            'connectivity': list(connectivity)}


if __name__ == '__main__':

    from timeit import default_timer as timer

    start = timer()
    data1 = reading_from_arc_file('QT4C/qt4c_dimer.txyz')
    data2 = reading_from_arc_file('QT4C/qt4c_dimer.arc')
    end = timer()
    print(end-start)

    start = timer()
    data3 = reading_from_arc_file_mmap('QT4C/qt4c_dimer.txyz')
    data4 = reading_from_arc_file_mmap('QT4C/qt4c_dimer.arc')
    end = timer()
    print(end-start)

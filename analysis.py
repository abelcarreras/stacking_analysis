#!/usr/bin/env python
__test__ = False

import descriptor
import iofile
import numpy as np
import json
import argparse

parser = argparse.ArgumentParser(description='stacking analysis options')

parser.add_argument('input_file', metavar='input_file', type=str,
                    help='JSON format input file')

parser.add_argument('molecule_file', metavar='trajectory', type=str,
                    help='Tinker XYZ or ARC file')

parser.add_argument('--relax', metavar='N', type=int, default=0,
                    help='relax steps [default: 0]')

parser.add_argument('--max', metavar='N', type=int, default=20000,
                    help='max steps [default: 20000]')

parser.add_argument('--cutoff', metavar='F', type=float, default=5.0,
                    help='cutoff pair interaction in Angstrom [default: 5.0]')

parser.add_argument('--sampling', metavar='N', type=int, default=10,
                    help='trajectory sampling [default: 10]')

args = parser.parse_args()


def average_angles(angle_list, symmetry=360):
    """
    Proper average of angles
    :param angle_list: list of angles to be averaged
    :param symmetry: set the symmetry
    :return: average
    """

    average = 0
    for i, angle in enumerate(angle_list):

        #print('*', angle, average)
        div = (angle - average + symmetry // 2) // symmetry
        angle -= div * symmetry
        #print(angle)

        if i == 0:
            average += angle
        else:
            average = (average + angle/float(i)) / (float(i+1)/i)

        div = average // symmetry
        average -= div * symmetry

    return average


def get_angle_between_vectors(v1, v2, symmetry=360, reflexion=False, n1=None, n2=None):
    """
    calculate angle between to vectors v1 & v2 in degrees
    """

    if n1 is not None and n2 is not None:
        pj1 = v1 - np.dot(v1, n2)*n2
        pj2 = v2 - np.dot(v2, n1)*n1
        angle = average_angles([np.arccos(np.dot(pj1, v2)) * 180 / np.pi,
                                np.arccos(np.dot(pj2, v1)) * 180 / np.pi],
                               symmetry=180)

    else:
        angle = np.arccos(np.dot(v1, v2)) * 180 / np.pi

    div = angle // symmetry
    angle -= div * symmetry

    if reflexion:
        mid_point = symmetry // 2
        if angle > mid_point:
            angle = symmetry - angle

    return angle


def get_distance_between_planes(center0, center1, normal0, normal1):
    """
    calculate distance between two planes defined by their centers and normal vectors
    """
    return np.average(np.abs([np.dot(normal0, center1 - center0),
                              np.dot(normal1, center0 - center1)]))


def get_sliding_between_planes(center0, center1, normal0, normal1, vector0, vector1, absolute=True):
    """
    calculate sliding between two planes defined by their centers and normal vectors
    """

    # displacement of planes in Cartesian coordinates
    disp_vector0 = center0 + normal0 * np.dot(normal0, center1 - center0) - center1
    disp_vector1 = center1 + normal1 * np.dot(normal1, center0 - center1) - center0

    # Orient in the same direction in order to be able to average
    disp_vector1 *= -1

    # Generate new basis to express the sliding displacement
    #    First basis vector as linear combination of monomer orientation vector
    ref_vectorb = (vector0 + vector1)/np.linalg.norm(vector0 + vector1)
    #    Second basis vector orthogonal to first basis vector and plane normal vectors
    vector0b = np.cross(ref_vectorb, normal0)
    vector1b = np.cross(ref_vectorb, normal1)

    #print(vector0b, vector1b)

    # Change basis to "parallel/perpendicular to new basis" coordinates
    disp_vector0b = np.array([np.dot(disp_vector0, ref_vectorb), np.dot(disp_vector0, vector0b)])
    disp_vector1b = np.array([np.dot(disp_vector1, ref_vectorb), np.dot(disp_vector1, vector1b)])

    #print('disp_vectorb(0 & 1):', disp_vector0b, disp_vector1b)

    #plt3d.quiver(*(center0).tolist() + (disp_vector0 * 15).tolist(), color='blue')
    #plt3d.quiver(*(center1).tolist() + (disp_vector1 * 15).tolist(), color='violet')

    #plt3d.quiver(*(center0).tolist() + (ref_vector * 15).tolist(), color='violet')
    #plt3d.quiver(*(center1).tolist() + (ref_vector * 15).tolist(), color='violet')

    #plt3d.quiver(*(center0).tolist() + (vector0b * 15).tolist(), color='green')
    #plt3d.quiver(*(center1).tolist() + (vector1b * 15).tolist(), color='green')

    # Average displacements from two monomers
    disp_vectorb = np.average([disp_vector0b, disp_vector1b], axis=0)

    # Get absolute value
    if absolute:
        disp_vectorb = np.abs(disp_vectorb)

    #print('disp_vectorb', disp_vectorb)
    return disp_vectorb


def get_significative_pairs(centers, normals, radius=5):

    import itertools

    num_monomers = len(centers)

    list_pairs = []
    for pair in itertools.combinations(range(num_monomers), 2):
        distance_planes = get_distance_between_planes(centers[pair[0]], centers[pair[1]],
                                               normals[pair[0]], normals[pair[1]])

        distance_centers = np.linalg.norm(centers[pair[1]] - centers[pair[0]])
        #print(pair, distance)
        if distance_planes < radius and distance_centers < radius:
            list_pairs.append(pair)
    return list_pairs


#data = iofile.reading_from_arc_file_mmap('QT4C/qt4c_helix.arc', nrelax=100, nmax=190)
data = iofile.reading_from_arc_file_mmap(args.molecule_file, nrelax=args.relax, nmax=args.max)


with open(args.input_file) as f:
    input_data = json.load(f)

monomer_atoms = input_data['monomer_atoms']
center_range = input_data['center_range']
align_atoms = input_data['align_atoms']

# create 3D plot
if __test__:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    plt3d = plt.figure().gca(projection='3d')
    plt3d.set_xlim3d(left=-20, right=20)
    plt3d.set_ylim3d(bottom=-15, top=20)
    plt3d.set_zlim3d(bottom=-5, top=20)

rotation_angle = []
distance = []
slides = []
neighbors = []
total_centers = []

for coordinates in data['trajectory'][::args.sampling]:
    total_atoms = len(coordinates)
    total_monomers = total_atoms/monomer_atoms

    vectors = []
    centers = []
    normals = []

    for i in range(total_monomers):
        coordinates_center = coordinates[center_range[0] + monomer_atoms*i:
                                         center_range[1] + monomer_atoms*i]
        center, normal, vector = descriptor.get_plane(coordinates_center, direction=align_atoms)
        # print(point, normal, vector)

        if __test__:
            #plt3d.scatter(*coordinates.T, color='orange')
            plt3d.scatter(*coordinates_center.T, color='blue', s=50)

            plt3d.scatter(*center, color='orange')
            plt3d.quiver(*center.tolist() + (vector * 15).tolist(), color='red')
            plt3d.quiver(*(center).tolist() + (normal * 5).tolist(), color='green')

        vectors.append(vector)
        centers.append(center)
        normals.append(normal)

    pairs = get_significative_pairs(centers, normals, radius=args.cutoff)
    neighbors.append(len(pairs))

    if len(pairs) == 0:
        continue

    # rotation_angle_i = []
    # distance_i = []
    # slides_i = []

    for pair in pairs:
        rotation_angle.append(get_angle_between_vectors(vectors[pair[0]], vectors[pair[1]],
                                                        symmetry=180, reflexion=True,
                                                        n1=normals[pair[0]], n2=normals[pair[1]]))  # in degrees

        distance.append(get_distance_between_planes(centers[pair[0]], centers[pair[1]],
                                                    normals[pair[0]], normals[pair[1]]))
        slides.append(get_sliding_between_planes(centers[pair[0]], centers[pair[1]],
                                                 normals[pair[0]], normals[pair[1]],
                                                 vectors[pair[0]], vectors[pair[1]]))

    # rotation_angle.append(average_angles(rotation_angle_i, symmetry=180))  # in degrees
    # distance.append(np.average(distance_i))
    # slides.append(np.average(slides_i, axis=0))
    total_centers.append(np.array(centers))

if __test__:
    plt.show()

np.savetxt('rotation.dat', rotation_angle)
np.savetxt('distance.dat', distance)
np.savetxt('slides.dat', slides)
np.savetxt('neighbors.dat', neighbors)

with file('total_centers.dat', 'w') as outfile:
    for slice_2d in np.array(total_centers):
        outfile.write('\n\n')
        np.savetxt(outfile, slice_2d)

# Plot data
#plt.title('Rotation angle')
#plt.ylim([0,90])
#plt.plot(rotation_angle, 'o')
#plt.show()

#plt.title('Rotation angle')
#plt.hist(rotation_angle, normed=1, bins=20)
#plt.show()

#plt.title('Distance')
#plt.ylim([0, 10])
#plt.plot(distance, 'o')
#plt.show()

#plt.title('Slide')
#plt.plot(slides, 'o')
#plt.show()

#plt.title('neighbors')
#plt.plot(neighbors)
#plt.show()

#plt.title('neighbors')
#plt.hist(neighbors, normed=1)
#plt.xticks(np.arange(-0.5, 8.5, 1.0), [str(k-1) for k in range(20)])
#plt.show()

#plt.hist(rotation_angle[500:], normed=1)



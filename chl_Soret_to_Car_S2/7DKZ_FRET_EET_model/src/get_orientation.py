
import numpy as np
import pandas as pd
import re
from collections import defaultdict

class AtomObj:
    def __init__(self, name, coord):
        self.name = name
        self.coord = np.array(coord)
        self.element = re.sub(r'[^A-Za-z]', '', name)[0]

class ResidueObj:
    def __init__(self, resname, resseq, chain):
        self.resname = resname
        self.resseq = resseq
        self.chain = chain
        self.atoms = {}

    def add_atom(self, atom_name, coord):
        self.atoms[atom_name] = AtomObj(atom_name, coord)

    def __getitem__(self, key):
        return self.atoms[key]

    def get_atoms(self):
        return self.atoms.values()
    
#-----------------------------------------------------
#-----------------------------------------------------

# Mass-weighted center of mass function
def get_center_of_mass(residue):
    masses = {'H':1.008, 'C':12.011, 'N':14.007, 'O':15.999, 'Mg':24.305, 'S':32.06}
    coords, weights = [], []
    for atom in residue.get_atoms():
        elem = atom.element.strip().title()
        m = masses.get(elem, 12.0)
        coords.append(atom.coord)
        weights.append(m)
    coords = np.array(coords)
    weights = np.array(weights)
    return np.average(coords, axis=0, weights=weights)

#-----------------------------------------------------
#-----------------------------------------------------


CHL_NAMES = {"CHL", "CLA", "CLB", "CHLA", "CHLB"}
CAR_NAMES = {"BCR", "LUT", "ZEA", "CAR", "PHY", "XAT"}

#--------------------------------------------------

def read_dat_structure(dat_file):

    residues = {}

    with open(dat_file) as f:
        for line in f:

            line = line.strip()
            if not line:
                continue

            m_res = re.search(r"Residue\s+(\w+)", line)
            m_seq = re.search(r"resseq=(\d+)", line)
            m_atom = re.search(r"<Atom\s+(\w+)>", line)
            m_coord = re.search(r"\[([^\]]+)\]", line)

            # skip lines that don't match expected pattern
            if not (m_res and m_seq and m_atom and m_coord):
                continue

            chain = line.split()[0]
            resname = m_res.group(1)
            resseq = int(m_seq.group(1))
            atom = m_atom.group(1)

            coord = np.array([float(x) for x in m_coord.group(1).split()])

            key = (chain, resname, resseq)

            if key not in residues:
                residues[key] = ResidueObj(resname, resseq, chain)

            residues[key].add_atom(atom, coord)

    return residues

#----------------------------------------------------


def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def kappa2(mu_d, mu_a, r_vec):
    r_hat = unit_vector(r_vec)
    return (np.dot(mu_d, mu_a) - 3 * np.dot(mu_d, r_hat) * np.dot(mu_a, r_hat))**2

# ------------------------------------------------------
# Chl-local axes and TDM vector
# ------------------------------------------------------
def define_chl_axes(residue):
    """
    Define local axes of a Chl residue:
    y-axis: along NB → ND
    x-axis: along NA → NC (perpendicular to y)
    z-axis: perpendicular to x & y
    """
    NB = residue["NB"].coord
    ND = residue["ND"].coord
    NA = residue["NA"].coord
    NC = residue["NC"].coord

    y_axis = ND - NB
    y_axis /= np.linalg.norm(y_axis)

    x_axis = NC - NA
    x_axis -= np.dot(x_axis, y_axis) * y_axis  # make perpendicular to y
    x_axis /= np.linalg.norm(x_axis)

    z_axis = np.cross(x_axis, y_axis)
    z_axis /= np.linalg.norm(z_axis)

    return x_axis, y_axis, z_axis

def tdm_vector_from_angle(angle_deg, x_axis, y_axis):
    """
    Rotate a unit vector along y-axis by angle (deg) in the x–y plane
    """
    theta = np.deg2rad(angle_deg)
    vec = np.cos(theta) * y_axis + np.sin(theta) * x_axis
    return vec / np.linalg.norm(vec)

# ------------------------------------------------------

def extract_chl_car_vectors(dat_file, angles_file, chl_tdm_label=None):

    angles_df = pd.read_csv(
        angles_file, sep="\t", header=None,
        names=["index", "TDMlabel", "type", "length", "angle"]
    )

    if chl_tdm_label is not None:
        filtered_df = angles_df[angles_df["TDMlabel"] == chl_tdm_label]
        if filtered_df.empty:
            raise ValueError(
                f"Chl TDM label '{chl_tdm_label}' not found."
            )
        angles_df = filtered_df

    residues = read_dat_structure(dat_file)

    chl_centers, car_centers = {}, {}
    chl_vectors_dict, car_vectors_dict = {}, {}

    for (chain, rname, resseq), residue in residues.items():

        res_id = f"{rname}{resseq}({chain})"

        # ----------------------------
        # Chlorophyll
        # ----------------------------
        if rname in CHL_NAMES:

            try:
                center = residue["MG"].coord
            except KeyError:
                center = get_center_of_mass(residue)

            chl_centers[res_id] = center

            x_axis, y_axis, z_axis = define_chl_axes(residue)

            row = angles_df.iloc[0]
            angle = row["angle"]

            vec = tdm_vector_from_angle(angle, x_axis, y_axis)

            chl_vectors_dict[(res_id, row["TDMlabel"])] = vec
            chl_vectors_dict[res_id] = [vec]

        # ----------------------------
        # Carotenoid
        # ----------------------------
        elif rname in CAR_NAMES:

            center = get_center_of_mass(residue)
            car_centers[res_id] = center

            try:
                start = residue["C6"].coord
                end = residue["C25"].coord
                vec = unit_vector(end - start)
                car_vectors_dict[res_id] = vec
            except KeyError:
                continue

    # convert to arrays
    chl_ids, chl_coords = zip(*chl_centers.items()) if chl_centers else ([], [])
    car_ids, car_coords = zip(*car_centers.items()) if car_centers else ([], [])

    return (
        np.array(chl_ids),
        np.array(chl_coords),
        chl_vectors_dict,
        np.array(car_ids),
        np.array(car_coords),
        car_vectors_dict
    )

#--------------------------------------------------------


# ------------------------------------------------------
# Compute kappa^2 matrix
# ------------------------------------------------------
def compute_kappa2_matrix(chl_ids, chl_coords, chl_vectors_dict,
                          car_ids, car_coords, car_vectors_dict):
    """
    Compute Chl x Car kappa^2 matrix
    """
    kappa2_matrix = np.zeros((len(chl_ids), len(car_ids)))

    for i, chl_id in enumerate(chl_ids):
        chl_vecs = chl_vectors_dict.get(chl_id, [])
        for j, car_id in enumerate(car_ids):
            car_vec = car_vectors_dict.get(car_id)
            if car_vec is None or len(chl_vecs) == 0:
                kappa2_matrix[i, j] = np.nan
                continue

            r_vec = car_coords[j] - chl_coords[i]
            k2_vals = [kappa2(mu_d, car_vec, r_vec) for mu_d in chl_vecs]
            kappa2_matrix[i, j] = np.mean(k2_vals)

    return kappa2_matrix

# ------------------------------------------------------



import numpy as np
import re
from collections import defaultdict


def extract_distance(target_coord_file):
    
    CHL_NAMES = {"CHL", "CLA", "CLB", "CHLA", "CHLB"}
    CAR_NAMES = {"BCR", "LUT", "ZEA", "CAR", "PHY", "XAT"}

    # atomic masses
    masses = {'H':1.008, 'C':12.011, 'N':14.007, 'O':15.999, 'Mg':24.305, 'S':32.06}

    # store residues
    residues = defaultdict(list)

    # ---- read file ----
    with open(target_coord_file) as f:
        for line in f:

            if not line.strip():
                continue

            chain = line.split()[0]

            m_res = re.search(r"Residue\s+(\w+)", line)
            m_seq = re.search(r"resseq=(\d+)", line)
            m_atom = re.search(r"<Atom\s+(\w+)>", line)
            m_coord = re.search(r"\[([^\]]+)\]", line)

            if not (m_res and m_seq and m_atom and m_coord):
                continue

            resname = m_res.group(1)
            resseq = int(m_seq.group(1))
            atom = m_atom.group(1)

            coord = np.array([float(x) for x in m_coord.group(1).split()])

            residues[(chain, resname, resseq)].append((atom, coord))


    # ---- output containers ----
    chl_res_ids = []
    chl_coords = []

    car_res_ids = []
    car_coords = []


    # ---- process residues ----
    for (chain, resname, resseq), atoms in residues.items():

        res_label = f"{resname}{resseq}({chain})"

        # ---------- Chlorophyll ----------
        if resname in CHL_NAMES:

            for atom, coord in atoms:
                if atom == "MG":
                    chl_res_ids.append(res_label)
                    chl_coords.append(coord)
                    break


        # ---------- Carotenoid ----------
        elif resname in CAR_NAMES:

            coords = []
            weights = []

            for atom, coord in atoms:

                elem = atom[0].upper()
                if atom.upper().startswith("MG"):
                    elem = "Mg"

                mass = masses.get(elem, 12.0)

                coords.append(coord)
                weights.append(mass)

            coords = np.array(coords)
            weights = np.array(weights)

            com = np.average(coords, axis=0, weights=weights)

            car_res_ids.append(res_label)
            car_coords.append(com)


    # ---- convert to arrays ----
    chl_coords = np.array(chl_coords)
    car_coords = np.array(car_coords)


    # ---- sorting helper ----
    def sort_key(label):
        m = re.match(r"(\w+)(\d+)\((\w)\)", label)
        resname, resseq, chain = m.groups()
        return (chain, int(resseq))


    # ---- sort chlorophylls ----
    chl_pairs = sorted(zip(chl_res_ids, chl_coords), key=lambda x: sort_key(x[0]))
    chl_res_ids, chl_coords = zip(*chl_pairs)

    chl_res_ids = np.array(chl_res_ids)
    chl_coords = np.array(chl_coords)


    # ---- sort carotenoids ----
    car_pairs = sorted(zip(car_res_ids, car_coords), key=lambda x: sort_key(x[0]))
    car_res_ids, car_coords = zip(*car_pairs)

    car_res_ids = np.array(car_res_ids)
    car_coords = np.array(car_coords)


    #print("Chlorophyll residues:", chl_res_ids)
    #print("Chlorophyll Mg coords:\n", chl_coords)
    print ('Num chl identified : ', len(chl_coords))
    #print("Carotenoid residues:", car_res_ids)
    #print("Carotenoid COM coords:\n", car_coords)
    print ('Num car identified : ', len(car_coords))
    
    return chl_res_ids, chl_coords, car_res_ids, car_coords
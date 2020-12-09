## Get header data lazily
Get chain and name.

    chains = []
    chain = {}
    for line in open('Mediator-final-v5.pdb'):
        if 'COMPND' not in line:
            pass
        elif 'MOL_ID' in line:
            if chain:
                chains.append(chain)
                chain = {}
        elif 'MODULE' in line:
            chain['module'] = re.search('MODULE\: (.*)\;', line).group(1)
        elif 'MOLECULE' in line:
            chain['name'] = re.search('MOLECULE\: (.*)\;', line).group(1)
        elif 'CHAIN' in line:
            chain['chain'] = re.search('CHAIN\: (.*)\;', line).group(1)
        else:
            raise ValueError(f'What is {line}')
    else:
        chains.append(chain)
    chains



## RMSD
This is the other way

    rmsd = pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric(ori)
    print(scorefxn_cart(mid), scorefxn_cart(ori), rmsd.calculate(mid))
    
    
## LocalRelax weights

These are the weights as seen in the [Wang et al 2016](https://elifesciences.org/articles/17219) paper

    # <Set scale_sc_dens_byres="E:0.56,D:0.56,R:0.76,K:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/>
    # is an utter mystery.
    weights = {"cart_bonded_length": 0.5,
               "cart_bonded_torsion": 0.5,
               "cart_bonded_angle": 1.0,
               "pro_close": 0.0,
               "fa_sol":0.0,
               "elec_dens_fast": 30,
               "rama": 0.0,
               "rama_prepro": 1.0}
    
    scorefxn_local = pyrosetta.create_score_function('ref2015_cart')
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    for name, value in weights.items():
        scorefxn_local.set_weight(stm.score_type_from_name(name), value)
        
This was not used in the final model.

## I-Tasser

Clean ITasser records.

    import os, re
    
    def clean_ter(folder):
        """
        I-Tasser adds stuff to TER records for complexes, PyMol does not like this.
        """
        for file in os.listdir(folder):
            if '.pdb' not in file:
                continue
            with open(os.path.join(folder, file)) as r:
                pdbblock=r.read()
            if not re.search('TER .*\n', pdbblock):
                continue
            with open(os.path.join(folder, file), 'w') as w:
                pdbblock=w.write(re.sub('TER .*\n', 'TER\n', pdbblock))

## Consurf

I don't recall which I used, so here are both...

Migrate Consurf scores from the file `consurf.grades` to a Pyrosetta pose.

    def add_bfactor_from_consurf(pose: pyrosetta.Pose, grades_filename:str):
        """
        Adds the bfactors from a consurf run to a pose based on the PDB residues.
        ``replace_res_remap_bfactors`` or ``set_bfactors_from_atom_id_map``
        were not used but may have been a cleaner strategy. This was quicker to write.
        """
        # parse file
        conscores = {} # key is MET1:A
        with open(grades_filename) as r:
            for row in r:
                row = row.strip()
                if not len(row) or not row[0].isdigit():
                    continue
                parts = row.split()
                conscores[parts[2]] = float(parts[3])
        # add to pose
        pdb_info = pose.pdb_info()
        pdb2pose = pdb_info.pdb2pose
        for con_res in conscores:
            pose_res = pdb2pose(res=int(con_res[3:-2]), chain=con_res[-1])
            assert pose.residue(pose_res).name3() == con_res[:3], f'{pose.residue(pose_res).name3()} ≠ {con_res[:3]}' 
            for i in range(pose.residue(pose_res).natoms()):
                pdb_info.bfactor(pose_res,i, conscores[con_res])
        return conscores
        
Or PDB file.

    def add_bfactor_from_consurf(pdb_filename: str, grades_filename:str):
        """
        Adds the bfactors from a consurf run to a pdb file based on the PDB residues.
        """
        # parse file
        conscores = {}  # key is MET1:A
        with open(grades_filename) as r:
            for row in r:
                row = row.strip()
                if not len(row) or not row[0].isdigit():
                    continue
                parts = row.split()
                conscores[parts[2]] = float(parts[3])
        # add to pdb file
        # pymol destroys metadata... so doing it manually
        # ATOM      1  N   MET A   1      85.911  68.442  35.731  1.00  0.00           N
        with open(pdb_filename) as r, open(pdb_filename.replace('.pdb', '.bfactor.pdb'), 'w') as w:
            for record in r:
                if record.find('ATOM') == 0:
                    resi = int(record[22:26].strip())
                    resn = record[17:20]
                    chain = record[21] if record[21] != ' ' else 'A'
                    bfactor = conscores[f'{record[17:20]}{resi}:{chain}']
                    # ATOM      1  N   MET A   1      85.911  68.442  35.731  1.00  0.00           N
                    #record = f'{record[:60]}{record[60:65]}{record[65:]}'
                    record = f'{record[:60]}{bfactor:> 5}{record[65:]}'
                w.write(record)
        return conscores
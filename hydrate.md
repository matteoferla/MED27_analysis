
## Hydrate

This took way way longer than it was meant to. Also hydration changes the PDB numbering which is an utter nuisance.

    # Touch first to make ehrm the magic work (?!)
    pyrosetta.rosetta.protocols.hydrate.Hydrate().apply(pose)
    
    pyrosetta.toolbox.generate_resfile.generate_resfile_from_pose(pose, 'pack.resfile', pack=False, design=False, input_sc=True, freeze=[], specific={})
    # check if works:
    # packtask = pyrosetta.standard_packer_task(pose)
    # pyrosetta.rosetta.core.pack.task.parse_resfile(pose, packtask, 'pack.resfile')
    
    pyrosetta.rosetta.protocols.hydrate.hydrate_hyfile(pose=pose,
                                                       hydrate_V=chain_vector,
                                                       resfile='pack.resfile')
                                                       
    # ================== action
    pyrosetta.rosetta.protocols.hydrate.Hydrate().apply(pose)
    
Waters were added as proven by:

    water_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
    water_selector.set_residue_name3('WAT')
    wv = water_selector.apply(pose)
    pyrosetta.rosetta.core.select.residue_selector.ResidueVector(wv)

The problem then is that the chain numbering has changed. So this is wrong:

    scorefxn = pyrosetta.create_score_function('ref2015')
    all_sele = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    constraint = pyrosetta.rosetta.protocols.constraint_generator.AtomPairConstraintGenerator()
    constraint.set_residue_selector(all_sele)
    constraint.set_ca_only(True)
    constraint.set_use_harmonic_function(True)
    constraint.set_max_distance(2.0)
    constraint.set_sd(0.5)
    add_csts = pyrosetta.rosetta.protocols.constraint_generator.AddConstraints()
    add_csts.add_generator(constraint)
    add_csts.apply(pose)
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
    
    chain_i = 21
    chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_i)
    pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    # chain_vector = chain_sele.apply(pose)
    and_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(chain_sele, water_selector)
    neigh_sele = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(selector=and_sele,
                                                                               distance=5,
                                                                               include_focus_in_subset=True)
    neigh_vector = neigh_sele.apply(pose)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(allow_bb=neigh_vector)
    movemap.set_chi(allow_chi=neigh_vector)
    movemap.set_jump(True)
    #I cannot seem to find a residue to upstream jump vector fxn
    # fold_tree.upstream_jump_residue(5) is broken.
    
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
    relax.set_movemap(movemap)
    relax.apply(pose)
    
    pose.dump_scored_pdb('mediator.f.l.med27.hydro2.pdb', scorefxn)




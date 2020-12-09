## Code snippets that are incorrect or took too long

## Per chain Cartesian relax

This took way way longer than it was meant to so was killed.

    scorefxn_cart = pyrosetta.create_score_function('ref2015_cart')
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
    scorefxn_cart.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
    
    for chain_i in range(1, pose.num_chains()+1):
        chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_i)
        chain_vector = chain_sele.apply(pose)
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=chain_vector)
        movemap.set_chi(allow_chi=chain_vector)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_cart, 1)
        relax.set_movemap(movemap)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.cartesian(True)
        relax.apply(pose)
        
    pose.dump_scored_pdb('mediator.fixed.cart_per_chain.pdb', scorefxn_cart)
    
## Hydrate

This took way way longer than it was meant to so was killed.



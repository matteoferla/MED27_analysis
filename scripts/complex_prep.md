## Start

    import pyrosetta
    from init_boilerplate import make_option_string
    pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                    ex1=None,
                                                    ex2=None,
                                                    mute='all',
                                                    ignore_unrecognized_res=True,
                                                    load_PDB_components=False,
                                                    ignore_waters=False)
                                   )
    scorefxn_cart = pyrosetta.create_score_function('ref2015_cart')
    scorefxn = pyrosetta.create_score_function('ref2015')
                                   
## Load

    ori = pyrosetta.Pose()
    filename = 'Mediator-final-v5.pdb'
    pyrosetta.rosetta.core.import_pose.pose_from_file(ori, filename)
    
    
## Fixed BB relax
A quick regular relax with fixed backbones to repack the sidechains —several added during import.
This does not use restrained CAs.

    movemap = pyrosetta.MoveMap()       
    scorefxn = pyrosetta.create_score_function('ref2015')
    movemap.set_bb(False)
    movemap.set_chi(True)
    movemap.set_jump(False)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
    relax.set_movemap(movemap)
    relax.apply(pose)
    print(scorefxn(pose))
    pose.dump_scored_pdb('mediator.fixed.pdb', scorefxn)
    
## LocalRelax

    #<LocalRelax name="local_rlx" scorefxn="dens" max_iter="100" ncyc="1" ramp_cart="0" K="16" nexp="2"/>
    relax = pyrosetta.rosetta.protocols.relax.LocalRelax()
    relax.set_sfxn(scorefxn_cart)
    relax.set_K(16)
    relax.set_max_iter(100)
    relax.set_ncyc(3)
    relax.set_nexp(2)
    print(scorefxn_cart(pose))
    relax.apply(pose)
    print(scorefxn_cart(pose))
    pose.dump_scored_pdb('mediator.fixed.local.pdb', scorefxn_cart)
    
## Per chain FastRelax

As the scoring will be done with the regular dihedral scorefunction, not a cartesian one, 
a restrained FastRelax per chain was done.

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
    
    for chain_i in range(1, pose.num_chains()+1):
        chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_i)
        chain_vector = chain_sele.apply(pose)
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=chain_vector)
        movemap.set_chi(allow_chi=chain_vector)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 1)
        relax.set_movemap(movemap)
        relax.apply(pose)
        
    pose.dump_scored_pdb('mediator.fixed.local2.per_chain.pdb', scorefxn)
    
## Humanise

The model is mouse, but I require 
## Humanise

The model is mouse, but I require human.


    Human>  MADVINVSVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREK
    Mouse>  MADVLSVGVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREK
    Diff.>  ....**.*..........................................
    
    Human>  AFIAHFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDK
    Mouse>  AFIANFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDK
    Diff.>  ....*.............................................
    
    Human>  TPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKA
    Mouse>  TPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKA
    Diff.>  ..................................................
    
    Human>  QPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKV
    Mouse>  QPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKV
    Diff.>  ..................................................
    
    Human>  IVVMRSLFIDRTIVKGYNENVYTEDGKLDIWSKSNYQVFQKVTDHATTAL
    Mouse>  IVVMRSLFIDRTIVKGYNESVYTEDGKLDIWSKSSYQVFQKVTDHATTAL
    Diff.>  ...................*..............*...............
    
    Human>  LHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFR
    Mouse>  LHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFR
    Diff.>  ..................................................
    
    Human>  TLEAFHDTCRQ
    Mouse>  TLEAFHDTCRQ
    Diff.>  ...........


This requires data generated in [metadata_assembly](metadata_assembly.md).

    import json
    
    with open('metadata.json', 'r') as r:
        metadata = json.load(r)

I believe there is some duplication in the code... but it works.

    from Bio import pairwise2
    from typing import Tuple
    
    def get_human2pose_alignment(entry: dict) -> Tuple[str, str]:
        """
        Returns human, mouse_pose
        """
        alignments = pairwise2.align.globalxs(entry['mouse_sequence'],
                                              entry['human_sequence'],
                                              -1, #open
                                              -0.1 #extend
                                             )
        al=alignments[0]
        # mouse to mouse.
        i = 0
        sequence = ''
        ref = iter(chain['pose_gap_sequence']+'-'*10_000)
        for resn in al[0]:
            if resn == '-':
                sequence += '-'
            else:
                sequence += next(ref)
        return (al[1], sequence) 
    
    def write_grishin(target_name, target_sequence, template_name, template_sequence, outfile):
        with open(outfile, 'w') as w:
            w.write(f'## {target_name} {template_name}\n')
            w.write(f'#\n')
            w.write('scores_from_program: 0\n')
            w.write(f'0 {target_sequence}\n')
            w.write(f'0 {template_sequence}\n')
            w.write('--\n')

The threading has to be done on each chain

    def get_chain_pose(pose, chain:int):
        """
        I suspect there is a simpler way.
        PyMOL might be safer.
        """
        copy = pose.clone()
        chain_end_res = pyrosetta.rosetta.core.pose.chain_end_res
        def chain_start_res(pose, chain:int):
            if chain == 1:
                return 1
            else:
                return pyrosetta.rosetta.core.pose.chain_end_res(pose, chain - 1) + 1
        keep = pyrosetta.rosetta.protocols.grafting.simple_movers.KeepRegionMover()
        # init accepts str. start/end accept str and int
        keep.start(chain_start_res(pose, chain))
        keep.end(chain_end_res(pose, chain))
        keep.apply(copy)
        return copy
        
The threading is done with the partial threading mover

    def thread(entry: dict, pose: pyrosetta.Pose) -> Tuple[pyrosetta.Pose, 
                                                           pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover]:
        human_seq, pose_seq = get_human2pose_alignment(entry)
        aln_file = f'{entry["human_uniprot_name"]}.aln'
        write_grishin(entry["human_uniprot_name"],
                      human_seq,
                      f'chain_{entry["chain"]}',
                      pose_seq,
                      aln_file)
        align = pyrosetta.rosetta.core.sequence.read_aln(format='grishin', filename=aln_file)
        chain_pose = get_chain_pose(pose, entry['number'])
        thread = pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover(align=align[1], 
                                                                        template_pose=chain_pose)
        query_pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.pose.make_pose_from_sequence(query_pose,
                                                            human_seq,
                                                           'fa_standard')
        thread.apply(query_pose)
        return query_pose, thread

This leaves the ring of unmoved residues.
I am not fully sure what the deal is with the "ring".
If the pose is copied beforehand the ring actually splits into smaller rings near each terminus/missing loop.
If I move the xyz coordinates it segfaults.

To fix this here are three fake selectors.

    
    class RingSelector:
        """
        Select all residues in the "ring"
        based upon 12A from origin
        NB> This is not actually a residue selector. The logical selectors will not accept it.
        """
        def __init__(self, radius=12):
            self.radius=radius
        
        def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
            sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
            for r in range(1, pose.total_residue() + 1):
                if abs(pose.residue(r).xyz(1).x) < 12:
                    sele[r] = 1
            return sele
        
    class AlteredSelector:
        """
        Select residues that were altered in the threading.
        NB> This is not actually a residue selector. The logical selectors will not accept it.
        """
        def __init__(self, threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover):
            self.threader = threader
        
        def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
            sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
            mapping = self.threader.get_qt_mapping(pose).mapping()
            for r in range(1, pose.total_residue() + 1):
                if mapping[r] == 0:
                    sele[r] = 1
            return sele
        
    class UnalteredSelector:
        """
        Select residues that were altered in the threading
        NB> This is not actually a residue selector. The logical selectors will not accept it.
        """
        def __init__(self, threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover):
            self.threader = threader
        
        def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
            sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
            mapping = self.threader.get_qt_mapping(pose).mapping()
            for r in range(1, pose.total_residue() + 1):
                if mapping[r] != 0:
                    sele[r] = 1
            return sele
            

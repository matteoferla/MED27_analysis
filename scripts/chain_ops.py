from Bio import pairwise2
from typing import Tuple, List, Union
import re


class ChainOps:
    """
    Works on the list of dict (metadata.json) with among others keys

    * number (pose number)
    * chain (chain letter)
    * gene_name (gene name)

    """

    def __init__(self, chains: List[dict]):
        self.chains = chains

    def get_entry_of_key(self, value, key):
        return [chain for chain in self.chains if chain[key] == value][0]

    def get_entry(self, value):
        """
        guess key...

        :param value:
        :return:
        """
        if isinstance(value, int):
            return self.get_entry_of_key(value, 'number')
        elif isinstance(value, str) and len(value) == 1:
            return self.get_entry_of_key(value, 'chain')
        elif isinstance(value, str):
            return self.get_entry_of_key(value, 'gene_name')
        elif isinstance(value, dict) and 'gene_name' in value:  # an entry was passed
            return value

    def get_pose_of_chain(self, pose, value):  # pose is a pyrosetta.Pose
        return pose.split_by_chain()[self.get_entry(value)['number']]

    def __getitem__(self, value):
        return self.get_entry(value)

    def __iter__(self):
        return self.chains


class Transmogrifier(ChainOps):
    """
    This is to convert a mutation from one species (``owned_label``) to another (``wanted_label``) based on provided sequences.

    """

    def __init__(self, chains: List[dict], wanted_label:str, owned_label: str):
        super().__init__(chains) #self.chains
        self.wanted_label = wanted_label
        self.owned_label = owned_label

    @classmethod
    def from_chain_ops(cls, chain_ops: ChainOps, wanted_label:str, owned_label: str):
        return cls(chain_ops.chains, wanted_label, owned_label)

    def align_seqs(self, chain_selection) -> Tuple[str, str]:
        """
        Align the human seq to the mouse and store it the chain dict
        """
        chain = self.chains[chain_selection]
        alignments = pairwise2.align.globalxs(chain[f'{self.wanted_label}_sequence'],
                                              chain[f'{self.owned_label}_sequence'],
                                              -1,  # open
                                              -0.1  # extend
                                              )
        al = alignments[0]
        chain[f'{self.wanted_label}_aln_sequence'] = al[0]
        chain[f'{self.owned_label}_aln_sequence'] = al[1]
        return al[0], al[1]

    def covert_A2B(self, seqA: str, seqB: str, resiA: int) -> int:
        """
        Given seqA and seqB as two gap aligned sequences,
        and an off-by-one residue number of seqA without counting gaps,
        return an off-by-one residue number of seqB without counting gaps.
        """
        assert resiA > 0, 'Negative number is no.'
        assert isinstance(resiA, int), 'Float is no'
        # first step
        get_resi2aln_map = lambda seq: [i for i, r in enumerate(seq) if r != '-']
        aln_pos = get_resi2aln_map(seqA)[resiA - 1]
        # second
        return get_resi2aln_map(seqB).index(aln_pos) + 1

    def transmogrify(self, mutation: str, chain_selection: Union[str, int, dict]) -> str:
        chain = self[chain_selection]
        human_resi = int(re.search(r'(\d+)', mutation).group(1))
        mouse_resi = self.covert_A2B(chain[f'{self.owned_label}_aln_sequence'],
                                     chain[f'{self.wanted_label}_aln_sequence'],
                                     human_resi)
        return re.sub(str(human_resi), str(mouse_resi), mutation)


class Murinizer(Transmogrifier):
    """
    Human --> Mouse
    """

    def __init__(self, chains: List[dict]):
        super().__init__(chains, 'mouse', 'human')


# from functools import partial
# med27_murinize = partial(murinize, entry=med27_chain)


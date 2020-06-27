from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

class MSAMap():

    def __init__(self, seqs: dict, gap_code="-"):
        self.seqs = seqs
        self.pos_msa_seq_map = {x: {} for x in self.seqs.keys()}
        self.pos_seq_msa_map = {x: {} for x in self.seqs.keys()}
        self.gap_code = gap_code

    def aln_len(self):
        return len(self.seqs[list(self.seqs)[0]])

    def samples(self):
        return list(self.seqs.keys())

    def init(self):
        curr_pos = {x: 0 for x in self.samples()}
        for msa_pos in range(self.aln_len()):
            for sample in self.samples():
                if self.seqs[sample][msa_pos] != self.gap_code:
                    self.pos_msa_seq_map[sample][msa_pos] = curr_pos[sample]
                    self.pos_seq_msa_map[sample][curr_pos[sample]] = msa_pos
                    curr_pos[sample] += 1

    def variants(self, refseq):
        variants = defaultdict(lambda: defaultdict(list))
        for msa_pos in range(self.aln_len()):
            if msa_pos in self.pos_msa_seq_map[refseq]:
                ref_pos = self.pos_msa_seq_map[refseq][msa_pos]
                ref_data = self.seqs[refseq][msa_pos]
                assert ref_data != "-"
                variant_id = f"{ref_data}_{ref_pos}"
                for sample in self.samples():
                    data = self.seqs[sample][msa_pos]
                    if ref_data == self.gap_code:
                        variants[variant_id]["del" + data].append(sample)
                    else:
                        variants[variant_id][data].append(sample)
            else:
                for sample in self.samples():
                    data = self.seqs[sample][msa_pos]
                    if data != self.gap_code:
                        variants[variant_id]["ins" + data].append(sample)



        return dict(variants)

    def pos_from_seq(self, seq_name_in, pos_in, seq_name_out):
        if pos_in in self.pos_seq_msa_map[seq_name_in]:
            msa_pos = self.pos_seq_msa_map[seq_name_in][pos_in]
            if msa_pos in self.pos_msa_seq_map[seq_name_out]:
                return self.pos_msa_seq_map[seq_name_out][msa_pos]
            else:
                raise ValueError(f"MSAPos {msa_pos} from {seq_name_out} not found")
        else:
            raise ValueError(f"Pos {pos_in} from {seq_name_in} not found")

    def subseq(self, seq_name_in, pos_in_start, pos_in_end, seq_name_out):
        seq = ""
        start = -1

        for pos_in in range(pos_in_start, pos_in_end):
            try:
                pos_out = self.pos_from_seq(seq_name_in, pos_in, seq_name_out)
                seq += self.seqs[seq_name_out][pos_out]
                if start == -1:
                    start = pos_out
            except ValueError:
                continue
        end = pos_out

        return SeqRecord(id=f'{seq_name_out}_{start}_{end}', name="", description="", seq=Seq(seq))

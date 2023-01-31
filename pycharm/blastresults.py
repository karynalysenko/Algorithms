from db import DB
from blast import SimpleBlast


class BlastResults:
    def __init__(self, query, w):
        self.query = query
        self.w = w
        self.sequences_obj = DB()
        self.ids_obj = DB()
        self.sequences_obj.extract_sequences()
        self.ids_obj.extract_ids()

    def process_results(self):
        # joins the best_hit to respective id_genebank (from database)
        all_seqs = []
        for num, seq in enumerate(self.sequences_obj.sequences):
            seq_id = []
            seq = ''.join(seq)
            blast = SimpleBlast(self.w, seq)
            results = blast.best_hit(self.query, seq, self.w)
            seq_id.append(results)
            seq_id.append(''.join(self.ids_obj.sequences[num]))
            all_seqs.append(seq_id)
        return all_seqs

    def sort_results(self):
        def sort_key(item):
            # sorting the tuple by the last number
            if item[0]:
                return item[0][-1]
            return 0

        return sorted(self.process_results(), key=sort_key, reverse=True)

    def __repr__(self):
        # organized printing
        print(
            '{:>9} {:>15} {:>18} {:>12}  {:>5}'.format('Query position', 'Seq position', 'Length of match', 'Max match',
                                                       'Id'))
        for i in self.sort_results():
            if i[0]:
                print('{:>6} {:>18} {:>15} {:>16}      |'.format(*i[0]), i[1])
        return '_______________________________________________________________'

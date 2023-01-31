class SimpleBlast:
    def __init__(self, w, seq):
        self.w = w
        self.seq = seq

    def query_map(self, query):
        # returns a dictionary where key is query of length w and values are
        # the indexes of query that match with key

        query_map = {}
        for i in range(len(query) - self.w + 1):
            window = query[i: i + self.w]
            if window in query_map:
                query_map[window].append(i)
            else:
                query_map[window] = [i]
        return query_map

    def hits(self, query_map, sequence):
        # returns list of tuples where 1st index is from query and 2nd is from seq
        hits = []
        for window_query, query_offsets in query_map.items():
            if window_query in sequence:
                seq_offset = sequence.find(window_query)
                while seq_offset != -1:  # because of find()
                    for query_offset in query_offsets:
                        hits.append((query_offset, seq_offset))
                    seq_offset = sequence.find(window_query, seq_offset + 1)  # because while loop and find()
        return hits

    def extend_hit_dir(self, query, sequence, q_offset, s_offset, direction):
        # direction: -1 or +1
        count = 0
        matches = 0
        while q_offset >= 0 and s_offset >= 0 and q_offset < len(query) and s_offset < len(sequence):
            count += 1
            if query[q_offset] == sequence[s_offset]:
                matches += 1
            if 2 * matches < count:
                break
            q_offset += direction
            s_offset += direction
        return q_offset - direction, s_offset - direction, matches, count

    def extend_hit(self, query, seq, hit, w):
        # given a hit, will extend the search of matches to left and right
        o1, o2 = hit
        left = self.extend_hit_dir(query, seq, o1 - 1, o2 - 1, -1)
        right = self.extend_hit_dir(query, seq, o1 + w, o2 + w, +1)

        O1, O2, ML, SL = left
        _, _, MR, SR = right

        #         print({'ML':ML, 'SR':SR, 'SL':SL,'MR':MR})
        return O1, O2, w + SL + SR, ML + w + MR

    def best_hit(self, query, seq, w):
        # finds best hit from all hits and respective extensions
        query_map = self.query_map(query)
        hits = self.hits(query_map, seq)
        best_hit = None
        best_score = 0
        for hit in hits:
            o1, o2, l, m = self.extend_hit(query, seq, hit, w)
            score = m - (l - m)
            if score > best_score or (score == best_score and l < best_hit[2]):
                best_hit = (o1, o2, l, m)
                best_score = score
        return best_hit
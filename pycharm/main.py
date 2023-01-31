if __name__ == '__main__':
    from blastresults import BlastResults

    query = "ACTG"
    w = 3
    blast_results = BlastResults(query, w)
    print(blast_results)



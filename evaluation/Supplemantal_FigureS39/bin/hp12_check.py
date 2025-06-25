def SimpleFastaParser(handle):
    #Migrated from from https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py can also process the afa file
    """Iterate over Fasta records as string tuples.
    Arguments:
     - handle - input stream opened in text mode
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip().upper())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def main(hp1_file, hp2_file):
    
    with open(hp1_file, "r") as f1, open(hp2_file, "r") as f2:
        f1_records = SimpleFastaParser(f1)
        f2_records = SimpleFastaParser(f2)

        total = 0
        same = 0

        for (t1, s1), (t2, s2) in zip(f1_records, f2_records):
            total += 1
            if s1 == s2:
                print(t1, t2)
                same += 1
        
        print(total, same, same/total)


main("chr1_HP1.fa", "chr1_HP2.fa")
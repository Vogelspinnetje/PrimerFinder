"""
PrimerFinder
8 februari 2024
Yesse Monkou

Dit programma zoekt alle mogelijke primers rond aangegeven gen op een
genoom. Vervolgens worden alle goede primers in een tekstbestand weergeven.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def genoom_inladen(bestands_naam: str, query: str) -> tuple[str, list]:
    sequence = ""
    gen_positie_ = []

    for record in SeqIO.parse(bestands_naam, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                if "locus_tag" in feature.qualifiers and feature.qualifiers[
                    "locus_tag"][0] == query:
                    sequence = record.seq
                    gen_positie_ = [feature.location.start,
                                    feature.location.end]

    return sequence, gen_positie_


def primers_maken(sequence: str, primer_lengtes: list[int], gc_ratio_: list[
    float], smeltpunt: list[int]) -> dict:
    primers = {"sequentie": [], "gc percentage": [], "smeltpunt": [],
                 "locatie": []}

    # Loop over de mogelijke lengtes van primers
    for lengte_primer in range(primer_lengtes[0], primer_lengtes[1]):
        for startpositie in range(len(sequence)):
            if startpositie > len(sequence) - lengte_primer:
                continue

            pos_primer = sequence[startpositie:startpositie + lengte_primer]
            if pos_primer[-1] not in "GC":
                continue

            aantal_nucleotiden = {nuc: pos_primer.count(nuc) for nuc in "ACGT"}
            gc_percentage = gc_fraction(Seq(pos_primer))
            if gc_ratio[0] > gc_percentage or gc_percentage > gc_ratio_[1]:
                continue

            berekend_smeltpunt = 2 * (aantal_nucleotiden['A'] +
                                    aantal_nucleotiden['T']) + 4 * (
                                    aantal_nucleotiden['G'] +
                                    aantal_nucleotiden['C'])
            if (smeltpunt[0] > berekend_smeltpunt or berekend_smeltpunt >=
                    smeltpunt[1]):
                continue

            primers["sequentie"].append(pos_primer)
            primers["gc percentage"].append(gc_percentage)
            primers["smeltpunt"].append(berekend_smeltpunt)

    return primers


def controleren_primer(primer, genbank_sequentie_, rev_com):
    for primer_sequentie in primer["sequentie"]:
        if rev_com:
            seq = primer_sequentie.reverse_complement()
        else:
            seq=primer_sequentie

        begin_pos = genbank_sequentie_.find(seq) + 1

        if begin_pos != -1:  # Als de sequentie gevonden is (find retourneert niet -1)
            # Bereken de eindpositie
            eind_pos = begin_pos + len(seq)
            primer["locatie"].append(f"{begin_pos} tot {eind_pos}")

    return primer


def print_primers(primers):
    tekst = []
    for primertjes in range(len(primers["sequentie"])):
        tekst.append(f"Sequentie: {primers['sequentie'][primertjes]}\n"
              f"GC percentage: {primers['gc percentage'][primertjes]}\n"
              f"Smeltpunt: {primers['smeltpunt'][primertjes]}\n"
              f"Locatie: {primers['locatie'][primertjes]}\n")

    return "\n".join(tekst)


if __name__ == "__main__":
    bestandsnaam = "sequence.gb"
    gezocht_gen = "lambdap07"
    afstanden = [100, 25]
    lengte = [17, 26]
    gc_ratio = [0.54, 0.60]
    smeltpunt_gebruiker = [55, 60]

    genbank_sequentie, gen_positie = genoom_inladen(bestandsnaam, gezocht_gen)
    seqf = genbank_sequentie[(gen_positie[0] - afstanden[0]):(gen_positie[0] -
                                                           afstanden[1])]
    seqb = Seq(genbank_sequentie[(gen_positie[1] + afstanden[1]):(
            gen_positie[1] + afstanden[0])])
    seqb = seqb.reverse_complement()

    f_primers = primers_maken(seqf, lengte, gc_ratio, smeltpunt_gebruiker)
    b_primers = primers_maken(seqb, lengte, gc_ratio, smeltpunt_gebruiker)
    f_primers = controleren_primer(f_primers, genbank_sequentie, False)
    b_primers = controleren_primer(b_primers, genbank_sequentie, True)

    print(f"Van gen {gezocht_gen}:")
    print(f"De forward primers:\n{print_primers(f_primers)}\n")
    print(f"De backward primers:\n{print_primers(b_primers)}")

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
import yaml


def load_config(config_pad: str) -> dict:
    with open(config_pad, "r") as config_yaml:
        return yaml.safe_load(config_yaml)


def genoom_inladen(bestands_naam: str, query: str) -> tuple[Seq, list[int, int]]:
    sequence: str = ""
    gen_positie_: list[int, int] = []

    for record in SeqIO.parse(bestands_naam, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "locus_tag" in feature.qualifiers and feature.qualifiers["locus_tag"][0] == query:
                    sequence: Seq = record.seq
                    gen_positie_: list[int, int] = [feature.location.start,
                                    feature.location.end]

    if sequence == "" or gen_positie_ == []:
        raise KeyError(f"Locus \"{query}\" is niet teruggevonden in "
                       f"opgegeven GenBank bestand. Zit deze locus "
                       f"daadwerkelijk in het bestand?")

    return sequence, gen_positie_


def primers_maken(sequence: Seq, primer_lengtes: list[int, int], gc_ratio_: list[float, float], smeltpunt: list[int, int]) -> dict:
    primers: dict = {"sequentie": [], "gc percentage": [], "smeltpunt": [], "locatie": []}

    for lengte_primer in range(primer_lengtes[0], primer_lengtes[1]):
        for startpositie in range(len(sequence)):
            if startpositie > len(sequence) - lengte_primer:
                continue

            pos_primer: Seq = sequence[startpositie:startpositie + lengte_primer]
            if pos_primer[-1] not in "GC":
                continue

            aantal_nucleotiden: dict = {nuc: pos_primer.count(nuc) for nuc in "ACGT"}
            gc_percentage: float = gc_fraction(Seq(pos_primer))
            if gc_ratio[0] > gc_percentage or gc_percentage > gc_ratio_[1]:
                continue

            berekend_smeltpunt: int = 2 * (aantal_nucleotiden['A'] +
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


def controleren_primer(primer: dict, genbank_sequentie_: str, rev_com: bool) -> dict:
    for primer_sequentie in primer["sequentie"]:
        if rev_com:
            seq: Seq = primer_sequentie.reverse_complement()
        else:
            seq: Seq = primer_sequentie

        begin_pos: int = genbank_sequentie_.find(seq) + 1

        if begin_pos != -1:
            eind_pos: int = begin_pos + len(seq)
            primer["locatie"].append(f"{begin_pos} tot {eind_pos}")

    return primer


def print_primers(primers: dict) -> str:
    tekst: list[str] = []
    
    if not primers["sequentie"]:
        return "-Geen primers gevonden-\n"

    for primertjes in range(len(primers["sequentie"])):
        tekst.append(f"Sequentie: {primers['sequentie'][primertjes]}\n"
              f"GC percentage: {primers['gc percentage'][primertjes]:.2f}\n"
              f"Smeltpunt: {primers['smeltpunt'][primertjes]}\n"
              f"Locatie: {primers['locatie'][primertjes]}\n")

    return "\n".join(tekst)


def main(afstanden: list[int, int], lengte: list[int, int], gc_ratio: list[float, float], smeltpunt: list[int, int], genbank_naam: str, gezocht_gen: str):
    genbank_sequentie: Seq
    gen_positie: list[int, int]
    genbank_sequentie, gen_positie = genoom_inladen(genbank_naam, gezocht_gen)
    
    seqf: Seq = genbank_sequentie[(gen_positie[0] - afstanden[0]):(gen_positie[0] - afstanden[1])]
    seqb: Seq = genbank_sequentie[(gen_positie[1] + afstanden[1]):(gen_positie[1] + afstanden[0])]
    seqb: Seq = seqb.reverse_complement()

    f_primers: dict = primers_maken(seqf, lengte, gc_ratio, smeltpunt_gebruiker)
    b_primers: dict = primers_maken(seqb, lengte, gc_ratio, smeltpunt_gebruiker)
    
    f_primers = controleren_primer(f_primers, genbank_sequentie, False)
    b_primers = controleren_primer(b_primers, genbank_sequentie, True)

    print(f"Van gen {gezocht_gen}:")
    print(f"De forward primers:\n{print_primers(f_primers)}\n")
    print(f"De backward primers:\n{print_primers(b_primers)}")
    

if __name__ == "__main__":
    config_pad: str = "config.yaml"
    genbank_naam: str = "sequence.gb"
    gezocht_gen: str = "lambdap07"

    configurations: dict = load_config(config_pad)
    afstanden: list[int, int] = configurations["primer_config"]["afstand"]
    lengte: list[int, int] = configurations["primer_config"]["lengte"]
    gc_ratio: list[float, float] = configurations["primer_config"]["gc_percentage"]
    smeltpunt_gebruiker: list[int, int] = configurations["primer_config"]["smeltpunt"]

    main(afstanden, lengte, gc_ratio, smeltpunt_gebruiker, genbank_naam, gezocht_gen)
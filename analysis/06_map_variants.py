"""Figures 6 and S12: conserved cysteine mapping and exact-allele evidence."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from _global_protein_alignment import alignment_primary, global_alignment, window13


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "inputs" / "variants"
REFERENCE = ROOT / "reference"


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reproduce_alignments(output):
    ledger = read_table(REFERENCE / "Fig6_mapping_sites.csv")
    # Keep the original age-response screen associated with the variant workflow.
    screened = [
        row for row in read_table(INPUT / "age_screen_input.tsv")
        if row["Mito_HighConf"] == row["Testable"] == "TRUE"
        and float(row["P.Value"]) < .05 and abs(float(row["logFC"])) > .5
    ]
    exact = [row for row in screened if row["N_CandidateSites"] == "1"]
    for row in exact:
        assert row["LocalizationStatus"] == "Single_site" and row["SiteAA"] == "C"
        assert row["CandidateSites"] == "C" + row["SiteLocation"]
    target_sites = {row["Accession"] + ":C:" + row["SiteLocation"] for row in exact}
    assert target_sites == {row["biological_key"] for row in ledger}
    assert len(exact) == 700 and len({(r["Omics"], r["FeatureID"]) for r in exact}) == 651
    groups = defaultdict(list)
    for row in ledger:
        if row["pair_key"] not in ("", "NA"):
            groups[row["pair_key"]].append(row)
    pairs = read_table(INPUT / "alignment_pairs.tsv")
    assert len(ledger) == 595 and sum(map(len, groups.values())) == 484
    assert len(pairs) == len(groups) == 313
    solved = []
    for index, pair in enumerate(pairs, 1):
        key = pair["pair_key"]
        mouse, human = pair["mouse_sequence"], pair["human_sequence"]
        rows = groups[key]
        centers = {int(row["mouse_position"]) for row in rows}
        result = global_alignment(mouse, human, sorted(centers), pair_cap=20_000_000)
        m_hash = hashlib.sha256(mouse.encode()).hexdigest()
        h_hash = hashlib.sha256(human.encode()).hexdigest()
        for row in rows:
            center = int(row["mouse_position"])
            assert mouse[center - 1] == "C" and row["mouse_sequence_sha256"] == m_hash
            assert row["human_sequence_sha256"] == h_hash
            positions = result["center_assignments"][str(center)]
            assert positions == json.loads(row["center_optimal_human_positions"]), row["biological_key"]
            assert result["score"] == float(row["alignment_score"]), row["biological_key"]
            assert all(abs(result[k] - float(row[k])) < 1e-8 for k in (
                "identity_percent", "mouse_coverage_percent", "human_coverage_percent"
            ))
            passed = alignment_primary(result["identical_pairs"], result["nongap_pairs"],
                                       len(mouse), len(human), positions)
            assert passed == (row["primary_mapping_pass"].lower() == "true")
            conserved = False
            if len(positions) == 1 and positions[0]:
                pos = positions[0]
                aa = human[pos - 1]
                assert int(row["human_position"]) == pos and row["human_AA_REF"] == aa
                assert row["wt_human_window"] == window13(human, pos)["wt_human_window"]
                conserved = passed and aa == "C"
            assert conserved == (row["primary_conserved_cys"].lower() == "true")
            solved.append({"biological_key": row["biological_key"], "pair_key": key,
                           "human_positions": json.dumps(positions), "primary_pass": passed,
                           "conserved_cysteine": conserved, "alignment_score": result["score"]})
        if index % 50 == 0:
            print(f"Verified {index}/{len(pairs)} alignments", flush=True)
    assert sum(row["conserved_cysteine"] for row in solved) == 444
    write_table(output / "conserved_cysteine_alignment.tsv", solved)


def reproduce_network(output):
    evidence = read_table(REFERENCE / "Fig6_network_exact_allele_witnesses.csv")
    edges = defaultdict(set)
    for row in evidence:
        assembly, chrom, pos, ref, alt = row["allele_key"].split(":")
        assert assembly == "GRCh38" and chrom and int(pos) > 0
        assert len(ref) == len(alt) == 1 and ref != alt
        edges[(row["event_key"], row["disease_id"])].add(row["allele_key"])
    expected = read_table(REFERENCE / "Fig6_network_event_disease_edges.csv")
    assert set(edges) == {(row["event_key"], row["disease_id"]) for row in expected}
    assert len(edges) == 821 and len({k[0] for k in edges}) == 594
    assert len({k[1] for k in edges}) == 168
    assert len({allele for alleles in edges.values() for allele in alleles}) == 591
    write_table(output / "variant_disease_network.tsv", [
        {"event_key": event, "disease_id": disease, "exact_alleles": json.dumps(sorted(alleles))}
        for (event, disease), alleles in sorted(edges.items())
    ])


def reproduce_flanking_contexts(output):
    relations = {r["relation_index"]: r for r in read_table(INPUT / "relation_coverage9190.tsv")}
    contexts = {r["context_key"]: r for r in read_table(INPUT / "original_candidate_contexts529.tsv")}
    projections = read_table(INPUT / "relation_context_projection73520.tsv")
    effect_rows = defaultdict(list)
    for row in read_table(INPUT / "mouse_corrected_effect_details.tsv"):
        effect_rows[(row["context_key"], row["background_effect_id"])].append(row)
    disease = defaultdict(set)
    for row in read_table(INPUT / "EXACT_ALLELE_DISEASE_EDGES.tsv"):
        disease[row["allele_key"]].add(row["disease_id"])
    selected, witnesses = [], set()
    for source in projections:
        if source["site_context_status"] != "NOMINAL_ELIGIBLE":
            continue
        index, key = source["relation_index"], source["context_key"]
        r = relations[index]
        if r["key"] not in disease:
            continue
        if r["signed_offset"] != r["corresponding_mouse_offset"] or r["aa_ref"] != r["mouse_residue_at_homologous_position"]:
            continue
        if not all(r[field] == "TRUE" for field in ("human_full13", "mouse_full13", "aligned_center_verified")):
            continue
        assert 0 < abs(int(r["signed_offset"])) <= 6
        assert r["molecular_consequence"] == "MISSENSE" and r["aa_ref"] != r["aa_alt"]
        effect_ids = json.loads(contexts[key]["nominal_eligible_effect_ids"])
        assert effect_ids
        for eid in effect_ids:
            effect, = effect_rows[(key, eid)]
            assert (float(effect["P.Value"]) < .05 and abs(float(effect["logFC"])) > .5 and
                    effect["source_link_kind"] == "EXACT_SOURCE" and effect["Contrast"] == "Old-Young")
        event = "|".join(r[c] for c in ("protein", "human_center_position", "aa_ref", "aa_position", "aa_alt"))
        selected.append({"relation_index": index, "context_key": key, "event_key": event,
                         "allele_key": r["key"], "site_key": r["biological_key"],
                         "signed_offset": r["signed_offset"]})
        for concept in disease[r["key"]]:
            witnesses.add((key, event, r["key"], concept))
    expected = read_table(REFERENCE / "Fig6_S12_flank_relation_contexts.csv")
    fields = ("relation_index", "context_key", "event_key", "allele_key", "site_key", "signed_offset")
    keys = lambda rows: {tuple(row[k] for k in fields) for row in rows}
    assert len(selected) == len(expected) == 603 and keys(selected) == keys(expected)
    assert len({(r["site_key"], r["event_key"]) for r in selected}) == 508
    assert len({r["site_key"] for r in selected}) == 136
    assert len({r["context_key"] for r in selected}) == 168
    assert len(witnesses) == 827
    write_table(output / "flanking_variant_contexts.tsv", selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "06_variants")
    args = parser.parse_args()
    output = args.output.resolve()
    for protected in (INPUT.parent, REFERENCE, ROOT / "analysis"):
        if output.is_relative_to(protected.resolve()):
            parser.error("choose an output folder outside analysis, inputs and reference")
    output.mkdir(parents=True, exist_ok=True)
    reproduce_network(output)
    reproduce_flanking_contexts(output)
    reproduce_alignments(output)
    print("Variant reproduction verified: 313 alignments, 444 conserved sites, 821 links, 508 flanking events.")


if __name__ == "__main__":
    main()

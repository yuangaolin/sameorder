#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Human–Gorilla 1:1 SD-block order analysis

Input files (current directory by default):
  result.txt
  ElementSDs_100.fasta
  blocks.fasta
  Genome_size_Indexes.txt
  BG_MosaicSDs.fasta      (optional but recommended)
  Clusters_100.fasta      (optional but recommended)

ElementSD coordinates use midStart / midEnd only.
"""

import csv
import os
import re
from collections import defaultdict

# -----------------------------
# inputs
# -----------------------------
RESULT_11 = "result.txt"
ELEMENT_FILE = "/home/fanghuixuan/SDCode/gor_joint/ElementSDs_100.fasta"
BLOCK_FILE = "/home/fanghuixuan/SDCode/gor_joint/blocks.fasta"
GENOME_INDEX = "/home/fanghuixuan/SDCode/gor_joint/Genome_size_Indexes.txt"
MOSAIC_FILE = "/home/fanghuixuan/SDCode/gor_joint/BG_MosaicSDs.fasta"
CLUSTER_FILE = "/home/fanghuixuan/SDCode/gor_joint/Clusters_100.fasta"

# -----------------------------
# outputs
# -----------------------------
OUT_PAIR_TABLE = "01_human_gor_1to1_block_pairs.tsv"
OUT_HUMAN_ORDER = "02_human_order.tsv"
OUT_GOR_ORDER = "03_gorilla_order.tsv"
OUT_CANDIDATES = "04_same_order_adjacent_pairs.tsv"
OUT_RUNS = "05_same_order_runs.tsv"
OUT_SUMMARY = "06_order_summary.txt"

HUMAN = "human"
GORILLA = "gorilla"


def get_species(header):
    low = header.lower()
    if "homo sapiens" in low:
        return HUMAN
    if "gorilla" in low:
        return GORILLA
    return "other"


def get_chr_name(header):
    m = re.search(r"chromosome\s+([0-9]+|X|Y)\b", header, re.I)
    if m:
        x = m.group(1).upper()
        return str(int(x)) if x.isdigit() else x
    if "mitochond" in header.lower():
        return "MT"
    return None


def chr_sort_key(x):
    if x is None:
        return (9, "")
    x = str(x)
    if x.isdigit():
        return (0, int(x))
    if x == "X":
        return (1, 23)
    if x == "Y":
        return (1, 24)
    if x == "MT":
        return (2, 25)
    return (3, x)


def read_genome_index(path):
    info = {}
    idx = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            info[idx] = {
                "species": get_species(line),
                "chr": get_chr_name(line),
                "header": line.strip(),
            }
            idx += 1
    return info


def read_elements(path, genome_info):
    """
    ElementSDs_100.fasta:
      index chr midStart midEnd midLength startC endC lengthC

    blocks.fasta uses signed 1-based element IDs, therefore eid=index+1.
    """
    elements = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.split()
            if len(p) < 8:
                continue
            try:
                index0 = int(p[0])
                joint_chr = int(p[1])
                mid_start = int(p[2])
                mid_end = int(p[3])
            except ValueError:
                continue

            eid = index0 + 1
            gi = genome_info.get(joint_chr, {})
            elements[eid] = {
                "eid": eid,
                "joint_chr": joint_chr,
                "species": gi.get("species"),
                "chr": gi.get("chr"),
                "start": mid_start,
                "end": mid_end,
            }
    return elements


def parse_block_ids(raw):
    return [int(x) for x in re.findall(r"[+-]?\d+", raw)]


def canonical_block(raw):
    return ",".join(str(x) for x in parse_block_ids(raw))


def read_1to1_result(path, elements):
    blocks = []
    seen = set()

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            raw = row["Block_Raw_Data"].strip()
            key = canonical_block(raw)
            if key in seen:
                continue
            seen.add(key)

            human = []
            gorilla = []

            for sid in parse_block_ids(raw):
                eid = abs(sid)
                x = elements.get(eid)
                if x is None:
                    continue
                rec = dict(x)
                rec["sign"] = "+" if sid > 0 else "-"
                if rec["species"] == HUMAN:
                    human.append(rec)
                elif rec["species"] == GORILLA:
                    gorilla.append(rec)

            if len(human) != 1 or len(gorilla) != 1:
                print(
                    "WARNING not 1 human : 1 gorilla:", raw,
                    "human=", len(human), "gorilla=", len(gorilla)
                )
                continue

            h = human[0]
            g = gorilla[0]
            blocks.append({
                "block_key": key,
                "raw": raw,
                "human_eid": h["eid"],
                "human_chr": h["chr"],
                "human_start": h["start"],
                "human_end": h["end"],
                "human_sign": h["sign"],
                "gorilla_eid": g["eid"],
                "gorilla_chr": g["chr"],
                "gorilla_start": g["start"],
                "gorilla_end": g["end"],
                "gorilla_sign": g["sign"],
                "orientation": "same" if h["sign"] == g["sign"] else "opposite",
            })

    return blocks


def assign_human_order(blocks):
    by_chr = defaultdict(list)
    for b in blocks:
        by_chr[b["human_chr"]].append(b)

    for chrom, arr in by_chr.items():
        arr.sort(key=lambda x: (x["human_start"], x["human_end"]))
        for rank, b in enumerate(arr, 1):
            b["human_rank"] = rank
            b["label"] = f"H{chrom}_{rank:06d}"
    return by_chr


def assign_gorilla_order(blocks):
    by_chr = defaultdict(list)
    for b in blocks:
        by_chr[b["gorilla_chr"]].append(b)

    for chrom, arr in by_chr.items():
        arr.sort(key=lambda x: (x["gorilla_start"], x["gorilla_end"]))
        for rank, b in enumerate(arr, 1):
            b["gorilla_rank"] = rank
    return by_chr


def classify_block(signed_ids, elements):
    human = []
    gorilla = []

    for sid in signed_ids:
        x = elements.get(abs(sid))
        if x is None:
            continue
        if x["species"] == HUMAN:
            human.append(x)
        elif x["species"] == GORILLA:
            gorilla.append(x)

    h = len(human)
    g = len(gorilla)

    if h == 1 and g == 1:
        typ = "1:1"
    elif h == 1 and g >= 2:
        typ = "1:N"
    elif h >= 2 and g == 1:
        typ = "N:1"
    elif h >= 2 and g >= 2:
        typ = "N:N"
    elif h >= 2 and g == 0:
        typ = "N:0"
    elif h == 1 and g == 0:
        typ = "1:0"
    elif h == 0 and g >= 2:
        typ = "0:N"
    elif h == 0 and g == 1:
        typ = "0:1"
    else:
        typ = "OTHER"

    return typ, human, gorilla


def read_all_blocks(path, elements):
    all_blocks = []
    assigned_eids = set()
    block_index = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            ids = parse_block_ids(line)
            if not ids:
                continue

            block_index += 1
            typ, human, gorilla = classify_block(ids, elements)
            assigned_eids.update(abs(x) for x in ids)
            all_blocks.append({
                "id": block_index,
                "raw": line,
                "type": typ,
                "human": human,
                "gorilla": gorilla,
            })

    return all_blocks, assigned_eids


def overlap(s1, e1, s2, e2):
    return max(s1, s2) < min(e1, e2)


def find_intervening_blocks(species, chrom, gap_start, gap_end, all_blocks, exclude_keys):
    hits = []
    if gap_end <= gap_start:
        return hits

    for b in all_blocks:
        if canonical_block(b["raw"]) in exclude_keys:
            continue

        records = b["human"] if species == HUMAN else b["gorilla"]
        for x in records:
            if x["chr"] != chrom:
                continue
            if overlap(gap_start, gap_end, x["start"], x["end"]):
                hits.append({
                    "block_id": b["id"],
                    "type": b["type"],
                    "eid": x["eid"],
                    "start": x["start"],
                    "end": x["end"],
                })
    return hits


def find_unassigned_elements(species, chrom, gap_start, gap_end, elements, assigned_eids):
    hits = []
    if gap_end <= gap_start:
        return hits

    for eid, x in elements.items():
        if eid in assigned_eids:
            continue
        if x["species"] != species or x["chr"] != chrom:
            continue
        if overlap(gap_start, gap_end, x["start"], x["end"]):
            hits.append(x)
    return hits


def read_mosaics(path, genome_info):
    mosaics = defaultdict(list)
    if not os.path.exists(path):
        print("WARNING missing:", path)
        return mosaics

    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 3:
                continue
            try:
                joint = int(p[0])
                start = int(p[1])
                end = int(p[2])
            except ValueError:
                continue

            gi = genome_info.get(joint, {})
            sp = gi.get("species")
            chrom = gi.get("chr")
            if sp in {HUMAN, GORILLA} and chrom:
                mosaics[(sp, chrom)].append((start, end))

    for key in mosaics:
        mosaics[key].sort()
    return mosaics


def find_mosaic(species, chrom, start, end, mosaics):
    for i, (s, e) in enumerate(mosaics.get((species, chrom), []), 1):
        if start >= s and end <= e:
            return i, s, e
    return None, None, None


def read_clusters(path, genome_info):
    clusters = defaultdict(list)
    if not os.path.exists(path):
        print("WARNING missing:", path)
        return clusters

    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.split()
            if len(p) < 7:
                continue
            try:
                cid = int(p[0])
                joint = int(p[1])
                start = int(p[2])
                end = int(p[3])
            except ValueError:
                continue

            gi = genome_info.get(joint, {})
            sp = gi.get("species")
            chrom = gi.get("chr")
            if sp in {HUMAN, GORILLA} and chrom:
                clusters[(sp, chrom)].append((cid, start, end))

    for key in clusters:
        clusters[key].sort(key=lambda x: x[1])
    return clusters


def find_clusters_between(species, chrom, start, end, clusters):
    if end < start:
        start, end = end, start

    hits = []
    for cid, s, e in clusters.get((species, chrom), []):
        if e < start:
            continue
        if s > end:
            break
        hits.append((cid, s, e))
    return hits


def block_hit_string(hits):
    if not hits:
        return "NONE"
    return ";".join(
        f"{x['type']}:E{x['eid']}:{x['start']}-{x['end']}"
        for x in hits
    )


def element_hit_string(hits):
    if not hits:
        return "NONE"
    return ";".join(
        f"E{x['eid']}:{x['start']}-{x['end']}"
        for x in hits
    )


def cluster_string(hits):
    if not hits:
        return "NONE"
    return ";".join(
        f"C{x[0]}:{x[1]}-{x[2]}"
        for x in hits
    )


def main():
    genome_info = read_genome_index(GENOME_INDEX)

    nh = sum(1 for x in genome_info.values() if x["species"] == HUMAN)
    ng = sum(1 for x in genome_info.values() if x["species"] == GORILLA)
    print("Genome index: human sequences =", nh, "gorilla sequences =", ng)

    if nh == 0:
        raise RuntimeError("Genome_size_Indexes.txt 中没有识别到 Homo sapiens")
    if ng == 0:
        raise RuntimeError("Genome_size_Indexes.txt 中没有识别到 Gorilla；请检查 header")

    elements = read_elements(ELEMENT_FILE, genome_info)
    blocks = read_1to1_result(RESULT_11, elements)
    human_order = assign_human_order(blocks)
    gorilla_order = assign_gorilla_order(blocks)
    all_blocks, assigned_eids = read_all_blocks(BLOCK_FILE, elements)
    mosaics = read_mosaics(MOSAIC_FILE, genome_info)
    clusters = read_clusters(CLUSTER_FILE, genome_info)

    # 1. basic 1:1 pair table
    with open(OUT_PAIR_TABLE, "w") as out:
        out.write(
            "Label\tBlock_Raw_Data\tHuman_Chr\tHuman_Start\tHuman_End\t"
            "Human_EID\tHuman_Sign\tHuman_Rank\t"
            "Gorilla_Chr\tGorilla_Start\tGorilla_End\tGorilla_EID\t"
            "Gorilla_Sign\tGorilla_Rank\tCross_Orientation\n"
        )
        for b in sorted(blocks, key=lambda x: (chr_sort_key(x["human_chr"]), x["human_start"])):
            out.write(
                f"{b['label']}\t{b['raw']}\t{b['human_chr']}\t{b['human_start']}\t"
                f"{b['human_end']}\t{b['human_eid']}\t{b['human_sign']}\t{b['human_rank']}\t"
                f"{b['gorilla_chr']}\t{b['gorilla_start']}\t{b['gorilla_end']}\t"
                f"{b['gorilla_eid']}\t{b['gorilla_sign']}\t{b['gorilla_rank']}\t"
                f"{b['orientation']}\n"
            )

    # 2. human order
    with open(OUT_HUMAN_ORDER, "w") as out:
        out.write(
            "Human_Chr\tHuman_Rank\tLabel\tHuman_Start\tHuman_End\t"
            "Gorilla_Chr\tGorilla_Start\tGorilla_End\tGorilla_Rank\n"
        )
        for chrom in sorted(human_order, key=chr_sort_key):
            for b in human_order[chrom]:
                out.write(
                    f"{chrom}\t{b['human_rank']}\t{b['label']}\t{b['human_start']}\t"
                    f"{b['human_end']}\t{b['gorilla_chr']}\t{b['gorilla_start']}\t"
                    f"{b['gorilla_end']}\t{b['gorilla_rank']}\n"
                )

    # 3. gorilla order
    with open(OUT_GOR_ORDER, "w") as out:
        out.write(
            "Gorilla_Chr\tGorilla_Rank\tLabel\tGorilla_Start\tGorilla_End\t"
            "Human_Chr\tHuman_Rank\tHuman_Start\tHuman_End\n"
        )
        for chrom in sorted(gorilla_order, key=chr_sort_key):
            for b in gorilla_order[chrom]:
                out.write(
                    f"{chrom}\t{b['gorilla_rank']}\t{b['label']}\t{b['gorilla_start']}\t"
                    f"{b['gorilla_end']}\t{b['human_chr']}\t{b['human_rank']}\t"
                    f"{b['human_start']}\t{b['human_end']}\n"
                )

    # 4. same-order adjacent pairs
    candidates = []
    total_human_edges = 0
    same_order_edges = 0
    reverse_order_edges = 0
    different_gorilla_chr = 0
    same_gorilla_chr_not_adjacent = 0

    for hchr, arr in human_order.items():
        for i in range(len(arr) - 1):
            a = arr[i]
            b = arr[i + 1]
            total_human_edges += 1

            if a["gorilla_chr"] != b["gorilla_chr"]:
                different_gorilla_chr += 1
                continue

            rank_diff = b["gorilla_rank"] - a["gorilla_rank"]
            if rank_diff == -1:
                reverse_order_edges += 1
            if rank_diff != 1:
                if rank_diff != -1:
                    same_gorilla_chr_not_adjacent += 1
                continue

            same_order_edges += 1

            human_gap_start = a["human_end"]
            human_gap_end = b["human_start"]
            gor_gap_start = a["gorilla_end"]
            gor_gap_end = b["gorilla_start"]

            exclude = {a["block_key"], b["block_key"]}

            human_other = find_intervening_blocks(
                HUMAN, a["human_chr"], human_gap_start, human_gap_end,
                all_blocks, exclude
            )
            gorilla_other = find_intervening_blocks(
                GORILLA, a["gorilla_chr"], gor_gap_start, gor_gap_end,
                all_blocks, exclude
            )

            human_unassigned = find_unassigned_elements(
                HUMAN, a["human_chr"], human_gap_start, human_gap_end,
                elements, assigned_eids
            )
            gorilla_unassigned = find_unassigned_elements(
                GORILLA, a["gorilla_chr"], gor_gap_start, gor_gap_end,
                elements, assigned_eids
            )

            hm1, _, _ = find_mosaic(
                HUMAN, a["human_chr"], a["human_start"], a["human_end"], mosaics
            )
            hm2, _, _ = find_mosaic(
                HUMAN, b["human_chr"], b["human_start"], b["human_end"], mosaics
            )
            gm1, _, _ = find_mosaic(
                GORILLA, a["gorilla_chr"], a["gorilla_start"], a["gorilla_end"], mosaics
            )
            gm2, _, _ = find_mosaic(
                GORILLA, b["gorilla_chr"], b["gorilla_start"], b["gorilla_end"], mosaics
            )

            same_human_mosaic = hm1 is not None and hm1 == hm2
            same_gorilla_mosaic = gm1 is not None and gm1 == gm2

            human_clusters = find_clusters_between(
                HUMAN, a["human_chr"], a["human_end"], b["human_start"], clusters
            )
            gorilla_clusters = find_clusters_between(
                GORILLA, a["gorilla_chr"], a["gorilla_end"], b["gorilla_start"], clusters
            )

            orientation_consistent = a["orientation"] == b["orientation"]

            human_other_types = {x["type"] for x in human_other if x["type"] != "1:1"}
            gorilla_other_types = {x["type"] for x in gorilla_other if x["type"] != "1:1"}
            other_sd_between = bool(human_other_types) or bool(gorilla_other_types)

            if other_sd_between:
                interpretation = "OTHER_SD_BETWEEN"
            elif human_unassigned or gorilla_unassigned:
                interpretation = "UNASSIGNED_ELEMENT_BETWEEN"
            elif not same_human_mosaic:
                interpretation = "DIFFERENT_HUMAN_MOSAIC"
            elif not same_gorilla_mosaic:
                interpretation = "DIFFERENT_GORILLA_MOSAIC"
            elif not orientation_consistent:
                interpretation = "ORIENTATION_INCONSISTENT"
            else:
                interpretation = "CLEAN_SAME_ORDER_CANDIDATE"

            candidates.append({
                "a": a,
                "b": b,
                "human_gap": human_gap_end - human_gap_start,
                "gorilla_gap": gor_gap_end - gor_gap_start,
                "human_other": human_other,
                "gorilla_other": gorilla_other,
                "human_unassigned": human_unassigned,
                "gorilla_unassigned": gorilla_unassigned,
                "human_clusters": human_clusters,
                "gorilla_clusters": gorilla_clusters,
                "same_human_mosaic": same_human_mosaic,
                "same_gorilla_mosaic": same_gorilla_mosaic,
                "orientation_consistent": orientation_consistent,
                "interpretation": interpretation,
            })

    # 5. candidate details
    with open(OUT_CANDIDATES, "w") as out:
        out.write(
            "Human_Chr\tGorilla_Chr\tLabel1\tLabel2\t"
            "Human_Rank1\tHuman_Rank2\tGorilla_Rank1\tGorilla_Rank2\t"
            "Human_EID1\tHuman_EID2\tGorilla_EID1\tGorilla_EID2\t"
            "Human_Coord1\tHuman_Coord2\tGorilla_Coord1\tGorilla_Coord2\t"
            "Human_Gap\tGorilla_Gap\tOrientation1\tOrientation2\t"
            "Orientation_Consistent\tSame_Human_Mosaic\tSame_Gorilla_Mosaic\t"
            "Human_Other_Blocks\tGorilla_Other_Blocks\t"
            "Human_Unassigned_Elements\tGorilla_Unassigned_Elements\t"
            "Human_Clusters_Between\tGorilla_Clusters_Between\tInterpretation\n"
        )

        for x in candidates:
            a = x["a"]
            b = x["b"]
            out.write(
                f"{a['human_chr']}\t{a['gorilla_chr']}\t{a['label']}\t{b['label']}\t"
                f"{a['human_rank']}\t{b['human_rank']}\t{a['gorilla_rank']}\t{b['gorilla_rank']}\t"
                f"{a['human_eid']}\t{b['human_eid']}\t{a['gorilla_eid']}\t{b['gorilla_eid']}\t"
                f"{a['human_start']}-{a['human_end']}\t{b['human_start']}-{b['human_end']}\t"
                f"{a['gorilla_start']}-{a['gorilla_end']}\t{b['gorilla_start']}-{b['gorilla_end']}\t"
                f"{x['human_gap']}\t{x['gorilla_gap']}\t{a['orientation']}\t{b['orientation']}\t"
                f"{x['orientation_consistent']}\t{x['same_human_mosaic']}\t{x['same_gorilla_mosaic']}\t"
                f"{block_hit_string(x['human_other'])}\t{block_hit_string(x['gorilla_other'])}\t"
                f"{element_hit_string(x['human_unassigned'])}\t{element_hit_string(x['gorilla_unassigned'])}\t"
                f"{cluster_string(x['human_clusters'])}\t{cluster_string(x['gorilla_clusters'])}\t"
                f"{x['interpretation']}\n"
            )

    # 6. maximal same-order runs
    candidate_edges = {(x["a"]["label"], x["b"]["label"]) for x in candidates}
    block_by_label = {b["label"]: b for b in blocks}
    runs = []

    for hchr, arr in human_order.items():
        current = []
        for i in range(len(arr) - 1):
            a = arr[i]
            b = arr[i + 1]
            edge = (a["label"], b["label"])

            if edge in candidate_edges:
                if not current:
                    current = [a["label"], b["label"]]
                else:
                    current.append(b["label"])
            else:
                if len(current) >= 2:
                    runs.append((hchr, current))
                current = []

        if len(current) >= 2:
            runs.append((hchr, current))

    with open(OUT_RUNS, "w") as out:
        out.write(
            "Human_Chr\tGorilla_Chr\tRun_Length\tLabels\t"
            "Human_Start\tHuman_End\tGorilla_Start\tGorilla_End\n"
        )
        for hchr, labels in runs:
            arr = [block_by_label[x] for x in labels]
            gor_chrs = {x["gorilla_chr"] for x in arr}
            gchr = next(iter(gor_chrs)) if len(gor_chrs) == 1 else "MULTI"
            out.write(
                f"{hchr}\t{gchr}\t{len(labels)}\t{','.join(labels)}\t"
                f"{arr[0]['human_start']}\t{arr[-1]['human_end']}\t"
                f"{arr[0]['gorilla_start']}\t{arr[-1]['gorilla_end']}\n"
            )

    # 7. summary
    interpretation_count = defaultdict(int)
    for x in candidates:
        interpretation_count[x["interpretation"]] += 1

    with open(OUT_SUMMARY, "w") as out:
        out.write(f"Total 1:1 human-gorilla blocks: {len(blocks)}\n")
        out.write(f"Human adjacent 1:1 pairs: {total_human_edges}\n")
        out.write(f"Different gorilla chromosome: {different_gorilla_chr}\n")
        out.write(f"Same gorilla chromosome but not adjacent: {same_gorilla_chr_not_adjacent}\n")
        out.write(f"Gorilla adjacent reverse-order pairs: {reverse_order_edges}\n")
        out.write(f"Same-order adjacent pairs: {same_order_edges}\n")
        out.write(f"Same-order runs: {len(runs)}\n")
        out.write("\nInterpretation counts:\n")
        for key in sorted(interpretation_count):
            out.write(f"{key}\t{interpretation_count[key]}\n")

    print("=" * 70)
    print("Total human-gorilla 1:1 blocks:", len(blocks))
    print("Human adjacent pairs:", total_human_edges)
    print("Same-order adjacent pairs:", same_order_edges)
    print("Reverse-order adjacent pairs:", reverse_order_edges)
    print("Different gorilla chromosome:", different_gorilla_chr)
    print("Same gorilla chr but not adjacent:", same_gorilla_chr_not_adjacent)
    print("Same-order runs:", len(runs))
    print("Interpretation counts:")
    for key in sorted(interpretation_count):
        print(" ", key, interpretation_count[key])
    print("=" * 70)
    print("Main output:", OUT_CANDIDATES)


if __name__ == "__main__":
    main()

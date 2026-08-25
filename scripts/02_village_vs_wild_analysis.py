#!/usr/bin/env python3
"""
Classify PASS variants in variants/village_wild_candidate_regions.PASS.vcf
as private/shared/core across the 11-sample pilot (5 village dogs, 6 wild
canids: 4 Grey Wolf, 1 Coyote; see metadata/pilot11_final.tsv), plus flag
sites where ALT carriage is fixed in exactly one of the two groups and
absent from the other -- candidate domestication-associated sites within
the 246 CDR set (Pendleton et al. 2018). n=5+6 pilot: hypothesis-
generating, not a population-genetics test.
"""
import csv
import subprocess
import sys
from pathlib import Path

VCF = Path("variants/village_wild_candidate_regions.PASS.vcf")
SAMPLE_SHEET = Path("metadata/pilot11_final.tsv")


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def carries_alt(gt):
    return any(a not in (".", "0") for a in gt.replace("|", "/").split("/"))


def main():
    pop_of = {}
    with open(SAMPLE_SHEET) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pop_of[row["Run"]] = row["Category"]

    samples = run(["bcftools", "query", "-l", str(VCF)]).strip().split("\n")
    gt_matrix = run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n", str(VCF)]
    ).strip().split("\n")

    village = [s for s in samples if pop_of.get(s) == "village"]
    wild = [s for s in samples if pop_of.get(s) == "wild"]

    private_count = {s: 0 for s in samples}
    shared_count = {s: 0 for s in samples}
    site_carrier_count = {}
    core_sites = []
    group_fixed_sites = []

    for line in gt_matrix:
        parts = line.split("\t")
        chrom, pos = parts[0], parts[1]
        gts = dict(zip(samples, parts[2:]))
        carriers = [s for s, gt in gts.items() if carries_alt(gt)]
        n = len(carriers)
        site_carrier_count[n] = site_carrier_count.get(n, 0) + 1
        if n == 1:
            private_count[carriers[0]] += 1
        elif n > 1:
            for s in carriers:
                shared_count[s] += 1
        if n == len(samples):
            core_sites.append((chrom, pos))

        vc = sum(1 for s in village if carries_alt(gts[s]))
        wc = sum(1 for s in wild if carries_alt(gts[s]))
        if (vc, wc) in ((len(village), 0), (0, len(wild))):
            group_fixed_sites.append((chrom, pos, vc, wc))

    total = sum(site_carrier_count.values())
    print(f"Total PASS sites: {total}  (samples: {len(samples)} = {len(village)} village + {len(wild)} wild)")
    print(f"  private (exactly 1 sample):  {site_carrier_count.get(1, 0)}")
    print(f"  shared  (2+ samples):        {total - site_carrier_count.get(0, 0) - site_carrier_count.get(1, 0)}")
    print(f"  core    (all {len(samples)} samples):    {len(core_sites)}")
    print()
    print(f"{'Sample':<14}{'Group':<10}{'Private':>10}{'Shared':>10}{'Total':>8}")
    for s in samples:
        print(f"{s:<14}{pop_of.get(s, '?'):<10}{private_count[s]:>10}{shared_count[s]:>10}{private_count[s] + shared_count[s]:>8}")

    print()
    print("Sites by number of carrier samples:")
    for n in sorted(site_carrier_count):
        print(f"  {n:>2} samples: {site_carrier_count[n]} sites")

    print()
    print(f"Group-fixed sites (ALT fixed in all village OR all wild, absent from the other, n={len(village)}+{len(wild)}): {len(group_fixed_sites)}")
    for chrom, pos, vc, wc in group_fixed_sites:
        which = "village-fixed" if vc == len(village) else "wild-fixed"
        print(f"  {chrom}:{pos}  {which}  (village={vc}/{len(village)}, wild={wc}/{len(wild)})")


if __name__ == "__main__":
    sys.exit(main())

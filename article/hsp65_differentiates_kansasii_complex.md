# A single hsp65 Sanger read differentiates the species of the *Mycobacterium kansasii* complex

*Draft, 2026-10-16. Numbers come from `output/hsp65_refalign/` (this repo) and from the reference-level analyses of mlsa-kansasii (`output/mlsa/main1_locus_discovery*`, `main4_hsp65_vs_wgs`). Authors and affiliations to be added.*

## Abstract

**Background.** The *M. kansasii* complex comprises at least seven species (*kansasii*, *persicum*, *pseudokansasii*, *innocens*, *attenuatum*, *ostraviense*, *gastri*) that differ in clinical relevance and drug susceptibility. Diagnostic laboratories commonly identify them by Sanger sequencing of a 441 bp hsp65 (groEL2) fragment. Recent genomic work has questioned whether hsp65, or any classic marker, is sufficient.

**Methods.** We tested hsp65 against whole-genome species (GTDB-Tk) in two ways. At reference level, we measured how many of the 21 species pairs hsp65 separates in 72 public genomes. In the clinical setting, we aligned every routine hsp65 Sanger read of 136 genome-typed isolates to one reference genome per species and called the species with the fewest differences, requiring at least 98% identity and at least two differences fewer than the next species. Three contaminated isolates were excluded.

**Results.** hsp65 separated 21 of 21 species pairs once three near-clonal genomes carrying an atypical hsp65 were set aside (within-species divergence up to 0.45%, between-species at least 1.36%). A perfect read (hsp65 extracted from each assembly) gave the genome species for 177 of 177 uncontaminated isolates. Of 133 analysed isolates, 109 received a call and **all 109 were correct** (197 of 197 individual reads; the 72 isolates with two or more callable reads never disagreed). The remaining 24 isolates got no call, and in none of them did a callable read point to the wrong species: 13 carried clean reads of *Rothia mucilaginosa*, 10 only unusable traces, and 1 only poor reads.

**Conclusion.** One good hsp65 read, interpreted by best match and margin rather than by a fixed identity cutoff, identifies the species of the *M. kansasii* complex. Failures are failures of the read, not of the marker, and are detectable by the same rules that make the correct calls.

## 1. Introduction

*M. kansasii* was long treated as one species with hsp65 PCR-restriction subtypes I to V [Picardeau 1997; Alcaide 1997]. Phylogenomics showed that the subtypes are separate species, with ANI between them of 88.4 to 94.2%, below the species cut-off [Tagini 2019], and *M. gastri* sits inside the group. Identification matters because the species differ in clinical relevance and resistance profile, and gastri and kansasii cannot be told apart by 16S, which is identical in them [Richter 1999]. Hsp65 sequencing on the 441 bp Telenti fragment [Telenti 1993] has been the routine alternative, with sequence databases and a 97% identity rule [McNabb 2004; Kim 2005]. Against this, a genomic study concluded that none of 16S, hsp65, rpoB, tuf and ITS distinguishes all species and recommended whole genomes [Jagielski 2019].

We argue that this conclusion is true for a fixed identity threshold and for the full set of public genomes taken at face value, but not for what a laboratory actually does with a read: align it, take the closest species and require a margin. We test that on genomes and on twelve years of routine reads (2014 to 2026).

## 2. Methods

**Truth.** Whole-genome species (GTDB-Tk, release 232 genomes as reference) of 180 clinical isolates sequenced on Illumina and assembled with Unicycler. CheckM contamination above the lab's flag marks an isolate as contaminated (Mkan329-035, -106, -113); these are excluded because a mixed sample has no single species.

**Reference level.** hsp65, 16S, ITS and six further loci were extracted from 72 GTDB genomes of the seven species. A pair of species is separated if the largest within-species distance is smaller than the smallest between-species distance (DNA-barcoding gap). We ran this with all genomes and without three Korean kansasii genomes (GCF_002705785.1, -825.1, -865.1) that carry a persicum-like hsp65 (see Discussion).

**Clinical reads.** All hsp65 `.ab1` traces in the lab archive (2014 to 2026; primers TB-11 and TB-12w) were quality-trimmed (modified Mott), assigned to their isolate through every TNR recorded for it, and aligned to one reference genome per species. Seeds of exact 15-mers locate the best window and a global alignment counts differences over the read. The call is the species with the fewest differences. A read is **ok** if its identity is at least 98% and the next-closest species has at least two more differences. An isolate is correct if it has at least one ok read and every ok read names the genome species. Reads without an ok call were also aligned to the isolate's own assembly; below 90% identity the read does not come from that genome. Two non-matching reads, one per sequence group, were identified by BLAST against NCBI core_nt.

**Ceiling.** To separate failures of the assay from limits of the marker, the hsp65 amplicon was also extracted in silico from each isolate's own assembly and classified with the same rules.

The code is `01_validate_hsp65_refalign.py` and `02_blast_reads.py`; the alignment module is vendored from mlsa-kansasii.

## 3. Results

### 3.1 At reference level hsp65 separates all species except through one lineage

| Locus | Species pairs separated (of 21) |
|---|---|
| gyrA, gyrB, rpoB, recA, tuf | 21 |
| **hsp65** | **17 (21 without three Korean genomes)** |
| ITS | 16 |
| 16S | 12 |
| secA1 | 6 |

All four failing hsp65 pairs involve the same three near-clonal genomes (kansasii with gastri, ostraviense, persicum, pseudokansasii). Without them, the largest within-species distance is 0.45% and the smallest between-species distance 1.36% over 441 bp, a gap of about 4 nucleotides. With them, within-kansasii divergence is 2.5%, which swallows the gap. 16S, in contrast, fails for reasons that no genome choice removes: several species have identical 16S sequences.

### 3.2 A perfect hsp65 read identifies every isolate

For all 177 uncontaminated isolates the in-silico amplicon had status ok and named the genome species (123 kansasii, 30 persicum, 20 pseudokansasii, 2 attenuatum, 1 innocens, 1 ostraviense). The amplicon was identical to its closest reference in 168 isolates, differed by one base in 4 and by two in 5. The closest other species was at least 4 differences further away (median 8).

### 3.3 Real Sanger reads: no call is common, a wrong call did not occur

Of 136 isolates with a Sanger hsp65 read, 133 remained after excluding contamination.

| Outcome | Isolates |
|---|---|
| Called, correct | **109** (kansasii 91, persicum 11, pseudokansasii 4, attenuatum 1, innocens 1, ostraviense 1) |
| Called, wrong | **0** |
| No call, reads are *Rothia mucilaginosa* | 13 |
| No call, only empty or very short traces | 10 |
| No call, only poor reads (93 to 94% identity) | 1 |

![Identity and margin of every hsp65 read per genome species](figures/hsp65_identity_margin.png)

**Figure 1.** Each point is one hsp65 read of an uncontaminated isolate, grouped by genome species. **A**, identity to the closest reference; filled points are `ok` reads (coloured by species), open grey points are reads without a call (44 reads: *Rothia*, poor and failed reads). **B**, number of differences between the closest and the next-closest species, for `ok` reads. Dashed lines are the calling thresholds. The two groups are separated by a gap: every `ok` read is at 98.26% or above, no other read exceeds 96.55%, and the smallest margin is 3 differences.

The 109 correct isolates had 197 ok reads (median 379 bp, 99.5% identity, median margin 8 differences, minimum 3). All 197 named the genome species. For the 72 isolates with two or more ok reads (forward and reverse, or repeats) the reads never disagreed.

### 3.4 What the 24 no-calls are

Thirteen isolates have clean traces that match no kansasii-complex reference (about 80% identity) and also not the isolate's own genome (75 to 80%). The reads form two sequence groups, one of 16 reads from 9 isolates sampled between 2015 and 2026 (94 to 100% identical), one shared by Mkan329-007 and -178. BLAST of one read per group gave ***Rothia mucilaginosa*, 98.9% and 97.6% identity over the whole read**. Rothia is an oral commensal, so a respiratory culture or specimen can plausibly carry it, and the TB-11/TB-12w primers evidently amplify its hsp65; we did not trace the source. The genomes of these isolates are clean (CheckM at most 3.3%) and their 16S and other reads match the complex species, so the culture is right and the hsp65 read is a contaminant read. The rules reject such reads without help: the identity is nowhere near 98%.

Ten isolates have only empty or sub-100 bp traces, and one (persicum) only reads at 93 to 94%. All 24 are cases for repeating the PCR.

## 4. Discussion

**Why hsp65 works here although "no marker distinguishes all species".** The statement of Jagielski et al. [2019] compares markers across species at the level of fixed thresholds and public genomes. We find two things that matter in practice. First, the marker fails through one identifiable lineage, not diffusely. Second, the call needs the best match and its margin, not an identity cut-off: kansasii and gastri are about 98.4% identical in hsp65, so a rule such as at least 97% [McNabb 2004] accepts both, whereas the closest species with a margin of at least two differences does not.

**The atypical lineage is real and must be flagged, not ignored.** The three Korean genomes are *M. kansasii* by ANI (98.5 to 99.3% to other kansasii, at most 94% to any other species) and identical to the type strain at gyrA, rpoB, 16S and other loci, yet carry a persicum-like hsp65 (11 SNPs in the amplicon; the "atypical hsp65 types" of Jagielski et al.). An hsp65 read of such a strain would be called persicum. Our rules keep these genomes as references, which lets a read that matches them be reported as kansasii with an atypical-hsp65 note and sent for a second locus. None of the clinical reads behaved like this. We therefore claim that hsp65 differentiates the complex, with this caveat made operational, and not that it is infallible.

**A second locus is a confirmation step, not a replacement.** gyrA, gyrB, rpoB, recA and tuf each separate all 21 pairs and would be better markers if one had to choose from scratch. The point is that hsp65 is already sequenced routinely, has a large reference base, and reaches the same answer on every isolate where a usable read exists. Adding gyrA is appropriate when the read hits the atypical lineage, the margin is below two, or the result must be certain for reporting.

**Mixed cultures and contaminant reads are visible.** The one discordant isolate in the data, Mkan329-035, is a mixed culture (kansasii and ostraviense, CheckM 113%); hsp65 correctly reported the dominant ostraviense component. Mixed templates can show as secondary peaks in the trace, and the *Rothia* reads were rejected by the identity rule.

**Limits.**
- The clinical set comes from one laboratory and is dominated by kansasii (105 of 133). The rarer species are represented by few isolates (pseudokansasii 13, persicum 12, innocens 1, attenuatum 1, ostraviense 1), and gastri has none; for it the claim rests on 2 reference genomes.
- The truth is the genome species. Isolates in which the genome and the Sanger sample could differ were handled only through the contamination flag and the own-genome check.
- Eighteen percent of isolates (24 of 133) had no usable call. A diagnostic claim has to include that the assay must be repeated for roughly one isolate in six.
- The rules were fixed on the same archive they are tested on, and the thresholds (98% identity, margin of two) were chosen on this data. A prospective batch is needed.
- Only two *Rothia* reads were BLASTed; the others were assigned by their identity to those.

## 5. Conclusion

For the *M. kansasii* complex, a single good hsp65 Sanger read, aligned to one reference genome per species and called by best match with a margin, gave the genome species for all 109 isolates that produced a callable read and for all 177 perfect in-silico amplicons. The marker fails only through a defined atypical lineage, which the same alignment can flag. The failures in practice were absent or contaminated reads, which the rules reject. Hsp65 sequencing is sufficient for routine species identification in this complex, with gyrA as the confirmation locus for the flagged cases.

## References

- Alcaide F, et al. J Clin Microbiol 1997;35:1959-1964. doi:10.1128/jcm.35.8.1959-1964.1997
- Jagielski T, et al. Genomic insights into the *Mycobacterium kansasii* complex: an update. Front Microbiol 2019;10:2918. doi:10.3389/fmicb.2019.02918
- Kim H, et al. Int J Syst Evol Microbiol 2005;55:1649-1656. doi:10.1099/ijs.0.63553-0
- McNabb A, et al. J Clin Microbiol 2004;42:3000-3011. doi:10.1128/JCM.42.7.3000-3011.2004
- Picardeau M, et al. J Clin Microbiol 1997;35:25-32. doi:10.1128/jcm.35.1.25-32.1997
- Richter E, et al. Identification of *Mycobacterium kansasii* by using a DNA probe (AccuProbe) and molecular techniques. J Clin Microbiol 1999;37:964-970. doi:10.1128/JCM.37.4.964-970.1999
- Tagini F, et al. Phylogenomics reveal that *Mycobacterium kansasii* subtypes are species-level lineages. Int J Syst Evol Microbiol 2019;69:1696-1704. doi:10.1099/ijsem.0.003378
- Telenti A, et al. Rapid identification of mycobacteria to the species level by PCR and restriction enzyme analysis. J Clin Microbiol 1993;31:175-178. doi:10.1128/jcm.31.2.175-178.1993

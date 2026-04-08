"""Pipeline module metadata and result parsing.

This module provides:
- Static metadata for each pipeline module (name, tool, description, etc.)
- Parsing of pipeline output files into structured results for the frontend.
"""

from __future__ import annotations

import csv
import json
import logging
import math
from pathlib import Path

logger = logging.getLogger("chatnav.pipeline")

# ── Static pipeline module definitions ──
# These match the actual pipeline steps the user sees.
# Numbering follows the README: Step 0, 8a, 8b, 9, 10, 11.

MODULES = [
    {
        "number": 1,
        "name": "Peptide Generation + MHC Binding",
        "subtitle": "From mutations to candidate peptides with binding predictions",
        "tool": "MHCflurry 2.2 + pyensembl",
        "description": (
            "Somatic mutations from the MAF/VCF are parsed and for each missense, "
            "frameshift, and in-frame indel, the real protein sequence is looked up "
            "from Ensembl GRCh38 (release 110). All 8-11mer peptide windows "
            "containing the mutation site are generated. MHCflurry "
            "Class1PresentationPredictor scores both mutant and wildtype peptides "
            "for binding affinity, antigen processing, and presentation against "
            "the patient's HLA alleles."
        ),
        "significance": (
            "This step converts raw mutation calls into actionable peptide "
            "candidates. Using real protein sequences (not translations of "
            "genomic coordinates) ensures accuracy. Running MHCflurry on "
            "wildtype peptides at the same positions enables agretopicity "
            "calculation — the differential binding that flags mutations "
            "creating novel immune targets."
        ),
        "outputs": [],
        "references": [
            "O'Donnell et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides. Cell Systems.",
        ],
        "output_dir": None,
    },
    {
        "number": 2,
        "name": "Immunogenicity + Foreignness Scoring",
        "subtitle": "Will the immune system actually react to this peptide?",
        "tool": "BigMHC-IM + k-mer proteome scan",
        "description": (
            "Three independent signals are computed for each peptide-HLA pair. "
            "BigMHC-IM is a transformer-based immunogenicity predictor (ensemble "
            "of 7 models) that scores whether a presented peptide will trigger T "
            "cell activation (0-1). Foreignness is measured by k-mer overlap "
            "against the full UniProt human proteome (20,659 proteins, 10.4M "
            "9-mers) — checking exact and 1-mismatch hits. Agretopicity is "
            "computed as log2(wt_binding_rank / mut_binding_rank): a high positive "
            "value means the mutation creates a new binding event."
        ),
        "significance": (
            "MHC binding is necessary but not sufficient — many bound peptides "
            "are simply ignored by T cells. BigMHC-IM predicts actual "
            "immunogenicity (AUC 0.92 on Gartner benchmark). Foreignness "
            "ensures peptides are sufficiently different from self-proteins "
            "to avoid immune tolerance."
        ),
        "outputs": ["immunogenicity_scores.tsv"],
        "references": [
            "Albert et al. (2023). BigMHC: accurate prediction of MHC-I antigen presentation and immunogenicity. Nature Machine Intelligence.",
        ],
        "output_dir": "module_08",
    },
    {
        "number": 3,
        "name": "Structural Scoring",
        "subtitle": "Is the mutated residue exposed to the T cell receptor?",
        "tool": "Position lookup / PANDORA / AlphaFold",
        "description": (
            "A tiered approach scores whether the mutated residue is exposed "
            "on the surface of the peptide-MHC complex where the T cell receptor "
            "can see it. Tier 1 (instant): position lookup table from crystal "
            "structure surveys — TCR-facing positions are P4-P7 for 9-mers. "
            "Tier 2 (minutes): PANDORA homology modelling generates a pMHC "
            "structure, then per-residue solvent-accessible surface area is "
            "computed at the mutation site. Tier 3 (hours): AlphaFold batch "
            "prediction via Tamarind.bio API for the top 5-10 candidates."
        ),
        "significance": (
            "A mutation buried inside the MHC groove cannot be recognised by "
            "T cells regardless of its immunogenicity score. Structural scoring "
            "penalises anchor-position mutations and rewards surface-exposed "
            "mutations at TCR-facing positions."
        ),
        "outputs": [],
        "references": [],
        "output_dir": None,
    },
    {
        "number": 4,
        "name": "Composite Ranking & Selection",
        "subtitle": "Choosing the best 20 neoantigens for the vaccine",
        "tool": "Weighted composite scorer + HLA diversity",
        "description": (
            "All scored candidates are ranked by a weighted composite formula: "
            "immunogenicity (0.35) + foreignness (0.15) + agretopicity (0.10) + "
            "binding (0.15) + stability (0.10) + expression (0.10) + CCF (0.05). "
            "Hard filters remove candidates with binding rank > 2%, TPM < 1.0, "
            "or CCF < 0.5. Bonuses are applied for frameshifts (+0.10), shared "
            "neoantigens (+0.15), and CD4 epitopes (+0.05). The top 20 are "
            "selected with HLA allele diversity enforcement."
        ),
        "significance": (
            "No single signal reliably predicts immunogenicity. The composite "
            "score balances multiple orthogonal signals to maximise the "
            "probability that at least one epitope triggers a clinically "
            "meaningful immune response. HLA diversity ensures the vaccine "
            "works across the patient's full HLA repertoire."
        ),
        "outputs": ["selected_neoantigens.tsv"],
        "references": [
            "Gartner et al. (2021). A machine learning model for ranking candidate HLA-I neoantigens. Nature Medicine.",
        ],
        "output_dir": "module_09",
    },
    {
        "number": 5,
        "name": "Polyepitope Design",
        "subtitle": "Assembling neoantigens into an optimised vaccine construct",
        "tool": "Greedy TSP + OptiVax linkers",
        "description": (
            "Selected neoantigens are assembled into a single polyepitope "
            "protein. A greedy travelling salesman algorithm optimises epitope "
            "ordering to minimise junctional neoepitope formation. AAY linkers "
            "separate MHC-I epitopes; GPGPG linkers separate MHC-II epitopes. "
            "A signal peptide is prepended for secretory pathway targeting, and "
            "the MITD domain is appended for efficient MHC-I loading."
        ),
        "significance": (
            "The ordering and linker design prevent the junctions between "
            "epitopes from creating unintended neoepitopes that could cause "
            "off-target immune responses. The signal peptide ensures the "
            "polyepitope enters the antigen processing pathway."
        ),
        "outputs": ["polyepitope.faa", "polyepitope_design.tsv"],
        "references": [],
        "output_dir": "module_10",
    },
    {
        "number": 6,
        "name": "mRNA Design",
        "subtitle": "Optimising codon usage and mRNA structure for translation",
        "tool": "LinearDesign + UTR library",
        "description": (
            "The polyepitope protein sequence is reverse-translated into "
            "optimised mRNA using LinearDesign for joint codon usage (CAI) and "
            "secondary structure (MFE) optimisation. The HBB 5' UTR and "
            "AES/mtRNR1 3' UTR are added, followed by a 120nt poly-A tail. "
            "All uridine positions are marked for N1-methylpseudouridine "
            "substitution, and CleanCap-AG co-transcriptional capping is "
            "specified."
        ),
        "significance": (
            "mRNA stability and translational efficiency directly affect "
            "how much protein each cell produces. LinearDesign finds the "
            "optimal tradeoff between codon adaptation and thermodynamic "
            "stability. Pseudouridine substitution reduces innate immune "
            "sensing of the mRNA itself."
        ),
        "outputs": ["mrna_sequence.fasta", "synthesis_spec.json"],
        "references": [
            "Zhang et al. (2023). LinearDesign: Efficient algorithms for optimized mRNA sequence design. Nature.",
        ],
        "output_dir": "module_11",
    },
]


def get_module_metadata() -> list[dict]:
    """Return a copy of the static module definitions."""
    return [dict(m) for m in MODULES]


def _determine_module_status(
    module_number: int, job_status: str, results_dir: Path
) -> dict:
    """Determine the status of a module based on its output files."""
    meta = next((m for m in MODULES if m["number"] == module_number), None)
    if meta is None:
        return {"status": "queued"}

    if job_status in ("pending",):
        return {"status": "queued"}

    # Check if module has an output directory with files
    output_dir = meta.get("output_dir")
    if output_dir:
        mod_dir = results_dir / output_dir
        has_outputs = mod_dir.exists() and any(mod_dir.iterdir()) if mod_dir.exists() else False
        if has_outputs:
            return {"status": "complete"}

    if job_status == "completed":
        return {"status": "complete"}
    if job_status == "running":
        return {"status": "queued"}
    return {"status": "queued"}


def build_modules_state(job_status: str, results_dir: Path) -> list[dict]:
    """Build the full module state list for a job."""
    modules = []
    for meta in MODULES:
        state = _determine_module_status(meta["number"], job_status, results_dir)
        modules.append({
            **meta,
            "status": state["status"],
            "started_at": None,
            "completed_at": None,
            "duration_seconds": None,
            "progress_pct": None,
            "progress_detail": None,
        })
    return modules


def parse_selected_neoantigens(tsv_path: Path) -> list[dict]:
    """Parse selected_neoantigens.tsv into a list of neoantigen dicts."""
    if not tsv_path.is_file():
        return []

    neoantigens = []
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader, 1):
            try:
                pep = row.get("peptide_sequence", row.get("peptide", ""))
                pep_len = len(pep)
                mut_pos = pep_len // 2
                tcr_positions = list(range(3, max(3, pep_len - 2)))

                neoantigens.append({
                    "rank": int(row["rank"]) if "rank" in row and row["rank"] else i,
                    "gene": row.get("gene", ""),
                    "mutation": row.get("mutation", row.get("aa_change", "")),
                    "peptide_sequence": pep,
                    "hla_allele": row.get("hla_allele", ""),
                    "mhc_class": row.get("mhc_class", "I"),
                    "composite_score": _float(row.get("composite_score", "0")),
                    "bigmhc_score": _float(row.get("bigmhc_immunogenicity", row.get("bigmhc_score", "0"))),
                    "foreignness_score": _float(row.get("foreignness", row.get("foreignness_score", "0"))),
                    "agretopicity": _float(row.get("agretopicity", "0")),
                    "binding_rank": _float(row.get("binding_rank", "0")),
                    "expression_tpm": _float(row.get("tpm", row.get("expression_tpm", "0"))),
                    "ccf": _float(row.get("ccf", "1.0")),
                    "structural_tier": row.get("structural_tier", "position"),
                    "structural_score": _float(row.get("structural_score", "0")),
                    "is_frameshift": row.get("is_frameshift", "").lower() in ("true", "1", "yes"),
                    "is_shared_neoantigen": row.get("is_shared", row.get("is_shared_neoantigen", "")).lower() in ("true", "1", "yes"),
                    "tcr_facing_positions": tcr_positions,
                })
            except Exception:
                logger.warning("Failed to parse neoantigen row %d", i, exc_info=True)
    return neoantigens


def parse_synthesis_spec(json_path: Path) -> dict | None:
    """Parse synthesis_spec.json for mRNA construct stats."""
    if not json_path.is_file():
        return None
    try:
        return json.loads(json_path.read_text())
    except Exception:
        logger.warning("Failed to parse synthesis spec", exc_info=True)
        return None


def build_mrna_segments(neoantigens: list[dict], spec: dict | None) -> list[dict]:
    """Build mRNA architecture segments for the map visualization."""
    n = len(neoantigens)
    if n == 0:
        return []

    total = spec.get("mrna_length_nt", 1200) if spec else 1200
    utr5_pct = 3.6
    utr3_pct = 9.5
    polya_pct = 9.6
    linker_pct = 0.5
    remaining = 100.0 - utr5_pct - utr3_pct - polya_pct - max(0, n - 1) * linker_pct
    epitope_pct = remaining / n if n > 0 else 0

    segments = [{"type": "utr5", "label": None, "length_nt": int(total * utr5_pct / 100), "pct_of_total": utr5_pct}]
    for i in range(n):
        segments.append({
            "type": "epitope",
            "label": str(i + 1),
            "length_nt": int(total * epitope_pct / 100),
            "pct_of_total": epitope_pct,
        })
        if i < n - 1:
            segments.append({
                "type": "linker",
                "label": None,
                "length_nt": int(total * linker_pct / 100),
                "pct_of_total": linker_pct,
            })
    segments.append({"type": "utr3", "label": None, "length_nt": int(total * utr3_pct / 100), "pct_of_total": utr3_pct})
    segments.append({"type": "polya", "label": None, "length_nt": int(total * polya_pct / 100), "pct_of_total": polya_pct})
    return segments


def _find_file(results_dir: Path, patient_id: str, *candidates: str) -> Path | None:
    """Try multiple paths to find a file, return first that exists."""
    for c in candidates:
        p = results_dir / c
        if p.is_file():
            return p
        p = results_dir / patient_id / c
        if p.is_file():
            return p
    return None


def _count_lines(path: Path) -> int:
    """Count non-header lines in a TSV."""
    with open(path) as f:
        return sum(1 for _ in f) - 1


def _count_unique_mutations(path: Path) -> int:
    """Count unique mutation IDs from immunogenicity_scores.tsv peptide_id column."""
    mutations: set[str] = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row.get("peptide_id", "")
            # peptide_id format: GENE_p.MUT_PEPTIDESEQ — strip last segment
            parts = pid.rsplit("_", 1)
            if len(parts) == 2:
                mutations.add(parts[0])
    return len(mutations)


# Gartner et al. 2021 benchmark calibration (Nature Medicine, n=61 NCI patients)
# Benchmark 3: 5 patients, 9 confirmed immunogenic in 100 top-20 slots
# PPV@20 = 9/100 = 0.09
# Recall@20 macro average = 0.37
BENCHMARK_PPV_PER_EPITOPE = 0.09
BENCHMARK_RECALL_AT_20 = 0.37
BENCHMARK_SOURCE = "Gartner et al. 2021 (Nature Medicine, n=61)"


def build_vaccine_result(results_dir: Path, patient_id: str) -> dict | None:
    """Build the full VaccineResult from pipeline output files."""
    # Find selected neoantigens
    neo_path = _find_file(results_dir, patient_id,
                          "module_09/selected_neoantigens.tsv")
    if neo_path is None:
        return None

    neoantigens = parse_selected_neoantigens(neo_path)
    if not neoantigens:
        return None

    # Parse synthesis spec
    spec_path = _find_file(results_dir, patient_id,
                           "module_11/synthesis_spec.json")
    spec = parse_synthesis_spec(spec_path) if spec_path else None

    n_neoantigens = len(neoantigens)

    # ── Construct stats from synthesis_spec.json (all REAL) ──
    mrna_length = spec.get("length_nt", spec.get("mrna_length_nt", 0)) if spec else 0
    gc_pct = spec.get("gc_content", spec.get("gc_pct", 0)) if spec else 0
    if gc_pct < 1:  # gc_content is 0-1 in the file
        gc_pct = gc_pct * 100
    cai = spec.get("lineardesign_cai", spec.get("cai", 0)) if spec else 0
    mfe = spec.get("lineardesign_mfe_kcal", spec.get("mfe_kcal_mol", 0)) if spec else 0
    polyepitope_aa = spec.get("protein_length_aa", 0) if spec else 0
    if polyepitope_aa == 0:
        polyepitope_aa = sum(len(n["peptide_sequence"]) for n in neoantigens)

    # ── Pipeline funnel stats (derive from real files where possible) ──
    # Try pipeline_stats.json first (written by rank_and_select.py)
    stats_path = _find_file(results_dir, patient_id, "pipeline_stats.json")
    pipeline_stats = {}
    if stats_path:
        try:
            pipeline_stats = json.loads(stats_path.read_text())
        except Exception:
            pass

    # n_windows: from pipeline_stats or count immunogenicity_scores.tsv rows
    imm_path = _find_file(results_dir, patient_id,
                          "module_08/immunogenicity_scores.tsv")
    n_windows = pipeline_stats.get("n_windows", 0)
    if n_windows == 0 and imm_path:
        n_windows = _count_lines(imm_path)

    # n_unique_mutations: from pipeline_stats or derive from peptide_ids
    n_variants = pipeline_stats.get("n_unique_mutations", 0)
    if n_variants == 0 and imm_path:
        n_variants = _count_unique_mutations(imm_path)

    # n_binders and n_filtered: from pipeline_stats if available
    n_binders = pipeline_stats.get("n_binders", 0)
    n_filtered = pipeline_stats.get("n_filtered", 0)

    segments = build_mrna_segments(neoantigens, spec)

    return {
        "n_neoantigens": n_neoantigens,
        "mrna_length_nt": mrna_length,
        "polyepitope_aa": polyepitope_aa,
        "gc_pct": gc_pct,
        "cai": cai,
        "mfe_kcal_mol": mfe,
        # Benchmark calibration — from Gartner et al. 2021 Benchmark 3
        # PPV@20 = 9 confirmed immunogenic / 100 top-20 slots = 0.09
        "ppv_per_epitope": BENCHMARK_PPV_PER_EPITOPE,
        "recall_at_20": BENCHMARK_RECALL_AT_20,
        "benchmark_source": BENCHMARK_SOURCE,
        # Pipeline funnel
        "n_variants": n_variants,
        "n_windows": n_windows,
        "n_binders": n_binders,
        "n_filtered": n_filtered,
        "neoantigens": neoantigens,
        "mrna_segments": segments,
    }


def _float(val: str) -> float:
    """Safe float conversion."""
    try:
        v = float(val)
        return 0.0 if math.isnan(v) else v
    except (ValueError, TypeError):
        return 0.0

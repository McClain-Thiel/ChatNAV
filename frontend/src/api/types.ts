export type JobStatus = 'pending' | 'running' | 'completed' | 'failed' | 'cancelled';
export type ModuleStatus = 'complete' | 'running' | 'queued' | 'failed';
export type MHCClass = 'I' | 'II';
export type StructuralTier = 'position' | 'pandora' | 'alphafold';
export type MRNASegmentType = 'utr5' | 'epitope' | 'linker' | 'utr3' | 'polya';

export interface JobSummary {
  id: string;
  patient_id: string;
  status: JobStatus;
  created_at: string;
  started_at: string | null;
  completed_at: string | null;
  error_message: string | null;
  params: Record<string, unknown>;
}

export interface ModuleState {
  number: number;
  name: string;
  subtitle: string;
  tool: string;
  status: ModuleStatus;
  started_at: string | null;
  completed_at: string | null;
  duration_seconds: number | null;
  progress_pct: number | null;
  progress_detail: string | null;
  outputs: string[];
  description: string;
  significance: string;
  references: string[];
  output_dir: string | null;
}

export interface JobDetail extends JobSummary {
  modules: ModuleState[];
  current_module: number;
  result: VaccineResult | null;
}

export interface Neoantigen {
  rank: number;
  gene: string;
  mutation: string;
  peptide_sequence: string;
  hla_allele: string;
  mhc_class: MHCClass;
  composite_score: number;
  bigmhc_score: number;
  foreignness_score: number;
  agretopicity: number;
  binding_rank: number;
  expression_tpm: number;
  ccf: number;
  structural_tier: StructuralTier;
  structural_score: number;
  is_frameshift: boolean;
  is_shared_neoantigen: boolean;
  tcr_facing_positions: number[];
}

export interface MRNASegment {
  type: MRNASegmentType;
  label: string | null;
  length_nt: number;
  pct_of_total: number;
}

export interface VaccineResult {
  n_neoantigens: number;
  mrna_length_nt: number;
  polyepitope_aa: number;
  gc_pct: number;
  cai: number;
  mfe_kcal_mol: number;
  ppv_per_epitope: number;
  recall_at_20: number;
  benchmark_source: string;
  n_variants: number;
  n_windows: number;
  n_binders: number;
  n_filtered: number;
  neoantigens: Neoantigen[];
  mrna_segments: MRNASegment[];
}

export interface OutputFile {
  path: string;
  size_bytes: number;
  modified_at: string;
}

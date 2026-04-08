import type { JobDetail } from '../api/types';
import { computeProbability } from '../api/helpers';
import SectionRule from './SectionRule';
import Tooltip from './Tooltip';
import MRNAMap from './MRNAMap';
import NeoantigenGrid from './NeoantigenGrid';
import PipelineProvenance from './PipelineProvenance';
import DownloadRow from './DownloadRow';

interface Props {
  job: JobDetail;
}

export default function ResultsDetail({ job }: Props) {
  const result = job.result;

  if (!result) {
    return (
      <div style={{
        maxWidth: 640,
        margin: '0 auto',
        padding: '80px 52px',
        textAlign: 'center',
      }}>
        <p style={{
          fontFamily: "'Spectral', serif",
          fontStyle: 'italic',
          fontSize: 16,
          color: 'var(--text-meta)',
        }}>
          Results will appear here when the pipeline completes.
        </p>
      </div>
    );
  }

  const probability = computeProbability(result.ppv_per_epitope, result.n_neoantigens);
  const probPct = Math.round(probability * 100);

  return (
    <div className="animate-reveal" style={{ maxWidth: 860, margin: '0 auto', padding: '44px 40px' }}>
      {/* Page header */}
      <div style={{ marginBottom: 32 }}>
        <div style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 9,
          letterSpacing: 2.5,
          textTransform: 'uppercase',
          color: 'var(--text-muted)',
          marginBottom: 4,
        }}>
          Vaccine Report
        </div>
        <h1 style={{
          fontFamily: "'EB Garamond', serif",
          fontSize: 32,
          fontWeight: 500,
          color: 'var(--text-primary)',
          marginBottom: 6,
        }}>
          {job.patient_id}
        </h1>
        <div style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 10,
          color: 'var(--text-meta)',
        }}>
          {(job.params.tumor_type as string) || ''} · {result.n_variants.toLocaleString()} variants
        </div>
        <div style={{ marginTop: 12 }}>
          <div style={{ borderTop: '2px solid var(--border-heavy)' }} />
          <div style={{ borderTop: '1px solid var(--border-mid)', marginTop: 4 }} />
        </div>
      </div>

      {/* Primary Outcome */}
      <SectionRule label="Primary Outcome" />
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 0, marginBottom: 32 }}>
        {/* Left: probability */}
        <div style={{
          borderLeft: '2px solid var(--border-heavy)',
          paddingLeft: 20,
          paddingRight: 16,
        }}>
          <Tooltip text="The estimated probability that at least one of the selected vaccine epitopes will trigger a T cell immune response. Based on the proportion of confirmed immunogenic peptides in the top 20 selections across 5 benchmark patients (Gartner et al. 2021). Does not account for immune escape, tumour microenvironment, or HLA loss.">
            <div className="t-label" style={{ marginBottom: 8, borderBottom: '1px dotted var(--text-muted)' }}>
              P(≥1 epitope immunogenic)
            </div>
          </Tooltip>
          <div style={{ display: 'flex', alignItems: 'baseline' }}>
            <span style={{
              fontFamily: "'EB Garamond', serif",
              fontSize: 52,
              fontWeight: 500,
              color: 'var(--text-primary)',
              lineHeight: 1,
            }}>
              {probPct}
            </span>
            <span style={{
              fontFamily: "'EB Garamond', serif",
              fontSize: 28,
              color: 'var(--text-meta)',
            }}>%</span>
          </div>
          <div style={{
            height: 4,
            background: 'var(--border-faint)',
            borderRadius: 2,
            marginTop: 8,
            marginBottom: 12,
          }}>
            <div style={{
              height: '100%',
              width: `${probPct}%`,
              background: 'var(--border-heavy)',
              borderRadius: 2,
              transition: 'width 0.8s ease',
            }} />
          </div>
          <div style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 12,
            color: 'var(--text-meta)',
            lineHeight: 1.5,
          }}>
            Simplified binomial model: P = 1 − (1 − p)<sup>n</sup> where p ≈ {result.ppv_per_epitope.toFixed(2)} per epitope, n = {result.n_neoantigens}. Assumes epitopes are independent. Does not account for immune escape, tumour microenvironment, HLA loss, or T cell exhaustion. Calibrated on Gartner et al. 2021 benchmark (PPV@20).
          </div>
        </div>

        {/* Right: benchmark */}
        <div style={{
          borderLeft: '1px solid var(--border-mid)',
          paddingLeft: 20,
        }}>
          <div className="t-label" style={{ marginBottom: 12 }}>
            Benchmark Calibration
          </div>
          {[
            { key: 'Source', value: result.benchmark_source, tip: 'The clinical validation dataset used to calibrate the probability estimate. 61 NCI patients with experimentally confirmed immunogenic and non-immunogenic neoantigens.' },
            { key: 'Recall@20', value: result.recall_at_20.toFixed(2), tip: 'Of all confirmed immunogenic peptides in the benchmark, 37% landed in the pipeline\'s top 20 selections. Higher is better — it means the ranking algorithm surfaces real hits.' },
            { key: 'PPV/epitope', value: result.ppv_per_epitope.toFixed(2), tip: 'Positive predictive value per epitope: the chance that any single selected epitope is truly immunogenic. Calculated as 9 confirmed positives across 100 top-20 slots in 5 benchmark patients.' },
            { key: 'Epitopes selected', value: String(result.n_neoantigens), tip: 'The number of neoantigen epitopes included in the vaccine construct. More epitopes increase the chance of at least one triggering a response, but each must meet quality thresholds.' },
          ].map(item => (
            <div key={item.key} style={{
              display: 'flex',
              justifyContent: 'space-between',
              padding: '6px 0',
              borderBottom: '1px solid var(--border-faint)',
              fontFamily: "'Inconsolata', monospace",
              fontSize: 9.5,
            }}>
              <Tooltip text={item.tip}>
                <span style={{ color: 'var(--text-muted)', borderBottom: '1px dotted var(--border-mid)' }}>{item.key}</span>
              </Tooltip>
              <span style={{ color: 'var(--text-secondary)' }}>{item.value}</span>
            </div>
          ))}
          <div style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 11,
            color: 'var(--text-faint)',
            marginTop: 8,
            lineHeight: 1.5,
          }}>
            CCF assumed 1.0 for MAF-only input. Expression integration partial.
          </div>
        </div>
      </div>

      {/* Vaccine Construct */}
      <SectionRule label="Vaccine Construct" />

      {/* Stats row */}
      <div style={{
        display: 'grid',
        gridTemplateColumns: 'repeat(6, 1fr)',
        borderTop: '1px solid var(--border-heavy)',
        borderBottom: '1px solid var(--border-mid)',
        marginBottom: 20,
      }}>
        {[
          { label: 'Neoantigens', value: result.n_neoantigens, tip: 'Number of neoantigen epitopes selected for the vaccine, chosen from all candidates after scoring and filtering.' },
          { label: 'mRNA Length', value: `${result.mrna_length_nt.toLocaleString()} nt`, tip: 'Total length of the synthesis-ready mRNA including 5\' UTR, coding sequence, 3\' UTR, and poly-A tail.' },
          { label: 'Polyepitope', value: `${result.polyepitope_aa} aa`, tip: 'Length of the encoded protein in amino acids. Includes all epitopes, linkers, signal peptide, and MITD domain.' },
          { label: 'GC Content', value: `${result.gc_pct.toFixed(1)}%`, tip: 'Percentage of guanine and cytosine bases in the mRNA. 50-60% is optimal for stability and translation in human cells.' },
          { label: 'CAI', value: result.cai.toFixed(3), tip: 'Codon Adaptation Index (0-1). Measures how well the codons match human tRNA abundance. Higher values mean faster translation. >0.8 is excellent.' },
          { label: 'MFE', value: `${result.mfe_kcal_mol.toFixed(1)}`, tip: 'Minimum Free Energy of the mRNA secondary structure in kcal/mol. More negative = more stable structure = longer half-life in the cell.' },
        ].map((stat, i) => (
          <div key={stat.label} style={{
            padding: '8px 12px',
            borderRight: i < 5 ? '1px solid var(--border-faint)' : 'none',
          }}>
            <Tooltip text={stat.tip}>
              <div style={{
                fontFamily: "'Inconsolata', monospace",
                fontSize: 8,
                letterSpacing: 2,
                textTransform: 'uppercase',
                color: 'var(--text-muted)',
                marginBottom: 4,
                borderBottom: '1px dotted var(--border-mid)',
                display: 'inline-block',
              }}>
                {stat.label}
              </div>
            </Tooltip>
            <div style={{
              fontFamily: "'EB Garamond', serif",
              fontSize: 16,
              fontWeight: 500,
              color: 'var(--text-primary)',
            }}>
              {stat.value}
            </div>
          </div>
        ))}
      </div>

      {/* mRNA Map */}
      <MRNAMap segments={result.mrna_segments} totalLength={result.mrna_length_nt} />

      {/* Selected Neoantigens */}
      <div style={{ marginTop: 32 }}>
        <SectionRule label="Selected Neoantigens" />
        <NeoantigenGrid neoantigens={result.neoantigens} />
      </div>

      {/* Pipeline Provenance */}
      <PipelineProvenance modules={job.modules} result={result} />

      {/* Download Row */}
      <DownloadRow jobId={job.id} />
    </div>
  );
}

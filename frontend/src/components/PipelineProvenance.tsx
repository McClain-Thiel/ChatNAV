import { useState } from 'react';
import type { ModuleState, VaccineResult } from '../api/types';
import { formatDuration } from '../api/helpers';

interface Props {
  modules: ModuleState[];
  result: VaccineResult;
}

export default function PipelineProvenance({ modules, result }: Props) {
  const [expanded, setExpanded] = useState(false);

  const funnelSteps = [
    { value: result.n_variants, label: 'mutations' },
    { value: result.n_windows, label: 'peptide windows' },
    ...(result.n_binders > 0 ? [{ value: result.n_binders, label: 'MHC binders' }] : []),
    ...(result.n_filtered > 0 ? [{ value: result.n_filtered, label: 'pass filters' }] : []),
    { value: result.n_neoantigens, label: 'selected' },
  ].filter(s => s.value > 0);

  return (
    <div style={{ marginTop: 32 }}>
      {/* Clickable header */}
      <button
        onClick={() => setExpanded(!expanded)}
        style={{
          display: 'flex',
          alignItems: 'center',
          width: '100%',
          gap: 8,
          background: 'none',
          border: 'none',
          padding: 0,
          cursor: 'pointer',
          marginBottom: expanded ? 20 : 0,
        }}
      >
        <div style={{ flex: 1, borderTop: '2px solid var(--border-heavy)' }} />
        <span style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 8.5,
          letterSpacing: 2.5,
          textTransform: 'uppercase',
          color: 'var(--text-meta)',
          whiteSpace: 'nowrap',
        }}>
          {expanded ? '▲' : '▼'} Pipeline Provenance
        </span>
        <div style={{ flex: 1, borderTop: '2px solid var(--border-heavy)' }} />
      </button>

      {expanded && (
        <div className="animate-reveal">
          {/* Funnel */}
          <div style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            gap: 6,
            marginBottom: 24,
            flexWrap: 'wrap',
          }}>
            {funnelSteps.map((step, i) => (
              <div key={step.label} style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
                {i > 0 && (
                  <span style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 14,
                    color: 'var(--border-mid)',
                  }}>→</span>
                )}
                <div style={{ textAlign: 'center' }}>
                  <div style={{
                    fontFamily: "'EB Garamond', serif",
                    fontSize: 20,
                    fontWeight: 600,
                    color: 'var(--text-primary)',
                  }}>
                    {step.value.toLocaleString()}
                  </div>
                  <div style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 8,
                    color: 'var(--text-muted)',
                    textTransform: 'uppercase',
                    letterSpacing: 1.5,
                  }}>
                    {step.label}
                  </div>
                </div>
              </div>
            ))}
          </div>

          {/* Module rows */}
          {modules.map(mod => (
            <div key={mod.number} style={{
              display: 'flex',
              alignItems: 'flex-start',
              gap: 12,
              padding: '8px 0',
              borderBottom: '1px solid var(--border-faint)',
            }}>
              {mod.status === 'complete' && mod.duration_seconds != null && (
                <span style={{
                  fontFamily: "'Inconsolata', monospace",
                  fontSize: 9,
                  color: 'var(--status-complete)',
                  background: 'var(--status-complete-bg)',
                  padding: '2px 6px',
                  borderRadius: 2,
                  whiteSpace: 'nowrap',
                }}>
                  ✓ {formatDuration(mod.duration_seconds)}
                </span>
              )}
              <div style={{ flex: 1, minWidth: 0 }}>
                <span style={{
                  fontFamily: "'EB Garamond', serif",
                  fontSize: 13.5,
                  fontWeight: 600,
                  color: 'var(--text-primary)',
                }}>
                  {mod.name}
                </span>
                <span style={{
                  fontFamily: "'Inconsolata', monospace",
                  fontSize: 9,
                  color: 'var(--text-muted)',
                  marginLeft: 8,
                }}>
                  {mod.tool}
                </span>
                {mod.outputs.length > 0 && (
                  <div style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 9,
                    color: 'var(--text-faint)',
                    marginTop: 2,
                  }}>
                    ↳ {mod.outputs.join(', ')}
                  </div>
                )}
              </div>
            </div>
          ))}

          {/* Footer */}
          <div style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 12,
            color: 'var(--text-faint)',
            marginTop: 16,
          }}>
            Pipeline version: ChatNAV v0.1.0 · Benchmark: {result.benchmark_source}
          </div>
        </div>
      )}
    </div>
  );
}

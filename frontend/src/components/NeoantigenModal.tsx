import { useState } from 'react';
import type { Neoantigen } from '../api/types';
import { getEpitopeColor, getScoreColor, getMutatedPosition, getTCRFacingPositions } from '../api/helpers';
import StructureViewer from './StructureViewer';

interface Props {
  neoantigen: Neoantigen;
  onClose: () => void;
}

const TIER_ICONS: Record<string, { icon: string; label: string; color: string }> = {
  alphafold: { icon: '◈', label: 'AlphaFold', color: 'var(--status-complete)' },
  pandora: { icon: '◆', label: 'PANDORA', color: '#4a3a1a' },
  position: { icon: '◇', label: 'Position lookup', color: 'var(--text-meta)' },
};

export default function NeoantigenModal({ neoantigen: neo, onClose }: Props) {
  const [showStructure, setShowStructure] = useState(false);
  const scoreColor = getScoreColor(neo.composite_score);
  const tier = TIER_ICONS[neo.structural_tier] ?? TIER_ICONS.position;
  const mutPos = getMutatedPosition(neo.peptide_sequence.length);
  const tcrPositions = neo.tcr_facing_positions.length > 0
    ? neo.tcr_facing_positions
    : getTCRFacingPositions(neo.peptide_sequence.length);

  return (
    <div
      style={styles.backdrop}
      onClick={onClose}
    >
      <div
        style={styles.modal}
        onClick={e => e.stopPropagation()}
      >
        {/* Close button */}
        <button onClick={onClose} style={styles.closeBtn}>×</button>

        {/* Header */}
        <div style={styles.header}>
          <div style={{
            width: 4,
            alignSelf: 'stretch',
            borderRadius: 2,
            background: getEpitopeColor(neo.rank),
            marginRight: 16,
            flexShrink: 0,
          }} />
          <div style={{ flex: 1 }}>
            <div style={{ display: 'flex', alignItems: 'baseline', gap: 8, marginBottom: 4 }}>
              <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 10, color: 'var(--text-faint)' }}>
                #{neo.rank}
              </span>
              <span style={{ fontFamily: "'EB Garamond', serif", fontSize: 22, fontWeight: 600, color: 'var(--text-primary)' }}>
                {neo.gene}
              </span>
              <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 12, color: 'var(--text-meta)' }}>
                {neo.mutation}
              </span>
              <div style={{ flex: 1 }} />
              {neo.is_frameshift && <Badge label="frameshift" color="var(--status-complete)" bg="var(--status-complete-bg)" />}
              {neo.is_shared_neoantigen && <Badge label="shared driver" color="#1a3a5c" bg="#eff6ff" />}
              <Badge label={`MHC Class ${neo.mhc_class}`} color="var(--text-meta)" bg="var(--bg-hover)" />
            </div>
            <div style={{ fontFamily: "'Inconsolata', monospace", fontSize: 10, color: 'var(--text-muted)' }}>
              {neo.hla_allele} · {tier.icon} {tier.label}
            </div>
          </div>
        </div>

        {/* Peptide sequence — large */}
        <div style={{ margin: '20px 0', textAlign: 'center' }}>
          <div style={{ fontFamily: "'Inconsolata', monospace", fontSize: 8, color: 'var(--text-muted)', letterSpacing: 2, textTransform: 'uppercase', marginBottom: 8 }}>
            Peptide Sequence
          </div>
          <div style={{ display: 'inline-flex', gap: 2 }}>
            {neo.peptide_sequence.split('').map((aa, i) => {
              const isMutated = i === mutPos;
              const isTCR = tcrPositions.includes(i);
              return (
                <span
                  key={i}
                  style={{
                    fontFamily: "'Inconsolata', monospace",
                    fontSize: 22,
                    fontWeight: isMutated ? 700 : 400,
                    color: isMutated ? '#8b1a1a' : isTCR ? 'var(--text-secondary)' : 'var(--text-muted)',
                    borderBottom: isMutated ? '2.5px solid #8b1a1a' : isTCR ? '1.5px solid var(--border-light)' : 'none',
                    padding: '0 2px',
                    lineHeight: 1.4,
                  }}
                >
                  {aa}
                </span>
              );
            })}
          </div>
          <div style={{ fontFamily: "'Inconsolata', monospace", fontSize: 8, color: 'var(--text-faint)', marginTop: 6 }}>
            <span style={{ color: '#8b1a1a', fontWeight: 700 }}>red</span> = mutated · <span style={{ color: 'var(--text-secondary)' }}>dark</span> = TCR-facing · <span style={{ color: 'var(--text-muted)' }}>light</span> = anchor
          </div>
        </div>

        <div style={{ borderTop: '2px solid var(--border-heavy)', marginBottom: 4 }} />
        <div style={{ borderTop: '1px solid var(--border-mid)', marginBottom: 20 }} />

        {/* Composite score bar */}
        <div style={{ marginBottom: 24 }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', marginBottom: 6 }}>
            <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 9, color: 'var(--text-muted)', letterSpacing: 1.5, textTransform: 'uppercase' }}>
              Composite Score
            </span>
            <span style={{ fontFamily: "'EB Garamond', serif", fontSize: 28, fontWeight: 500, color: scoreColor }}>
              {neo.composite_score.toFixed(3)}
            </span>
          </div>
          <div style={{ height: 4, background: 'var(--border-faint)', borderRadius: 2 }}>
            <div style={{
              height: '100%',
              width: `${neo.composite_score * 100}%`,
              background: scoreColor,
              borderRadius: 2,
              transition: 'width 0.5s ease',
            }} />
          </div>
        </div>

        {/* Score breakdown grid */}
        <div style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(4, 1fr)',
          gap: 0,
          border: '1px solid var(--border-light)',
          borderRadius: 3,
          overflow: 'hidden',
          marginBottom: 24,
        }}>
          {[
            { label: 'BigMHC-IM', value: neo.bigmhc_score, fmt: 3 },
            { label: 'Foreignness', value: neo.foreignness_score, fmt: 3 },
            { label: 'Agretopicity', value: neo.agretopicity, fmt: 2 },
            { label: 'Binding Rank', value: neo.binding_rank, fmt: 4 },
            { label: 'Expression', value: neo.expression_tpm, fmt: 1, suffix: ' TPM' },
            { label: 'CCF', value: neo.ccf, fmt: 2 },
            { label: 'Structural', value: neo.structural_score, fmt: 3 },
            { label: 'Tier', value: null, text: tier.label },
          ].map((item, i) => (
            <div key={item.label} style={{
              padding: '10px 12px',
              borderRight: (i % 4 !== 3) ? '1px solid var(--border-faint)' : 'none',
              borderBottom: i < 4 ? '1px solid var(--border-faint)' : 'none',
              background: i < 4 ? 'var(--bg-inset)' : 'transparent',
            }}>
              <div style={{
                fontFamily: "'Inconsolata', monospace",
                fontSize: 8,
                letterSpacing: 1.5,
                textTransform: 'uppercase',
                color: 'var(--text-muted)',
                marginBottom: 4,
              }}>
                {item.label}
              </div>
              <div style={{
                fontFamily: "'EB Garamond', serif",
                fontSize: 16,
                fontWeight: 500,
                color: 'var(--text-primary)',
              }}>
                {item.text ?? `${item.value!.toFixed(item.fmt)}${item.suffix ?? ''}`}
              </div>
            </div>
          ))}
        </div>

        {/* 3D Structure button + viewer */}
        <button
          onClick={() => setShowStructure(s => !s)}
          style={{
            display: 'block',
            width: '100%',
            padding: '10px 0',
            fontFamily: "'Inconsolata', monospace",
            fontSize: 10,
            letterSpacing: 1,
            color: showStructure ? 'var(--text-primary)' : 'var(--text-meta)',
            background: showStructure ? 'var(--bg-inset)' : 'transparent',
            border: '1px solid var(--border-mid)',
            borderRadius: 2,
            cursor: 'pointer',
            marginBottom: showStructure ? 0 : 16,
          }}
        >
          ◈ {showStructure ? 'Hide' : 'View'} pMHC Structure Prediction
        </button>

        {showStructure && (
          <div style={{ marginBottom: 16 }}>
            <StructureViewer
              peptide={neo.peptide_sequence}
              hlaAllele={neo.hla_allele}
              mutatedPosition={mutPos}
              tcrFacingPositions={tcrPositions}
            />
          </div>
        )}
      </div>
    </div>
  );
}

function Badge({ label, color, bg }: { label: string; color: string; bg: string }) {
  return (
    <span style={{
      fontFamily: "'Inconsolata', monospace",
      fontSize: 8.5,
      color,
      background: bg,
      border: `1px solid ${color}`,
      borderRadius: 2,
      padding: '2px 6px',
      lineHeight: 1,
    }}>
      {label}
    </span>
  );
}

const styles: Record<string, React.CSSProperties> = {
  backdrop: {
    position: 'fixed',
    inset: 0,
    background: 'rgba(28, 20, 16, 0.4)',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    zIndex: 1000,
  },
  modal: {
    background: 'var(--bg-card)',
    border: '1px solid var(--border-mid)',
    borderRadius: 4,
    padding: '28px 32px',
    maxWidth: 900,
    width: '70vw',
    maxHeight: '90vh',
    overflowY: 'auto',
    position: 'relative',
    boxShadow: '0 8px 32px rgba(28,20,16,0.15)',
    animation: 'reveal 0.15s ease',
  },
  closeBtn: {
    position: 'absolute',
    top: 12,
    right: 16,
    background: 'none',
    border: 'none',
    fontFamily: "'EB Garamond', serif",
    fontSize: 24,
    color: 'var(--text-muted)',
    cursor: 'pointer',
    lineHeight: 1,
    padding: 4,
  },
  header: {
    display: 'flex',
    marginBottom: 16,
  },
};

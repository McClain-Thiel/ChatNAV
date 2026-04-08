import type { Neoantigen } from '../api/types';
import { useHover } from '../context/HoverContext';
import { getEpitopeColor, getScoreColor, getMutatedPosition, getTCRFacingPositions } from '../api/helpers';

interface Props {
  neoantigen: Neoantigen;
  onSelect: (neo: Neoantigen) => void;
}

const TIER_ICONS: Record<string, { icon: string; color: string }> = {
  alphafold: { icon: '◈', color: 'var(--status-complete)' },
  pandora: { icon: '◆', color: '#4a3a1a' },
  position: { icon: '◇', color: 'var(--text-meta)' },
};

export default function NeoantigenCard({ neoantigen: neo, onSelect }: Props) {
  const { hoveredRank, setHoveredRank } = useHover();
  const isHighlighted = hoveredRank === neo.rank;
  const isDimmed = hoveredRank !== null && hoveredRank !== neo.rank;
  const scoreColor = getScoreColor(neo.composite_score);
  const tier = TIER_ICONS[neo.structural_tier] ?? TIER_ICONS.position;
  const mutPos = getMutatedPosition(neo.peptide_sequence.length);
  const tcrPositions = neo.tcr_facing_positions.length > 0
    ? neo.tcr_facing_positions
    : getTCRFacingPositions(neo.peptide_sequence.length);

  return (
    <div
      onClick={() => onSelect(neo)}
      onMouseEnter={() => setHoveredRank(neo.rank)}
      onMouseLeave={() => setHoveredRank(null)}
      style={{
        border: `1px solid ${isHighlighted ? 'var(--border-heavy)' : 'var(--border-light)'}`,
        borderLeft: `3px solid ${getEpitopeColor(neo.rank)}`,
        borderRadius: 3,
        padding: '10px 12px',
        background: isHighlighted ? 'var(--bg-hover)' : 'var(--bg-card)',
        opacity: isDimmed ? 0.6 : 1,
        cursor: 'pointer',
        transition: 'all 0.15s ease',
        boxShadow: isHighlighted ? '0 2px 8px rgba(28,20,16,0.08)' : 'none',
      }}
    >
      {/* Top row: rank, gene, mutation, badges */}
      <div style={{ display: 'flex', alignItems: 'baseline', gap: 6, marginBottom: 8 }}>
        <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 9, color: 'var(--text-faint)' }}>
          #{neo.rank}
        </span>
        <span style={{ fontFamily: "'EB Garamond', serif", fontSize: 14, fontWeight: 600, color: 'var(--text-primary)' }}>
          {neo.gene}
        </span>
        <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 9.5, color: 'var(--text-meta)' }}>
          {neo.mutation}
        </span>
        <div style={{ flex: 1 }} />
        <div style={{ display: 'flex', gap: 4 }}>
          {neo.is_frameshift && <SmallBadge label="fs" color="var(--status-complete)" />}
          {neo.is_shared_neoantigen && <SmallBadge label="shared" color="#1a3a5c" />}
          <SmallBadge label={`MHC-${neo.mhc_class}`} color="var(--text-meta)" />
        </div>
      </div>

      {/* Peptide sequence */}
      <div style={{ marginBottom: 8, letterSpacing: 1 }}>
        {neo.peptide_sequence.split('').map((aa, i) => {
          const isMutated = i === mutPos;
          const isTCR = tcrPositions.includes(i);
          return (
            <span
              key={i}
              style={{
                fontFamily: isMutated ? "'EB Garamond', serif" : "'Inconsolata', monospace",
                fontSize: 12,
                fontWeight: isMutated ? 700 : 400,
                color: isMutated ? '#8b1a1a' : isTCR ? 'var(--text-secondary)' : 'var(--text-muted)',
                borderBottom: isMutated ? '1.5px solid #8b1a1a' : isTCR ? '1px solid var(--border-light)' : 'none',
                padding: '0 0.5px',
              }}
            >
              {aa}
            </span>
          );
        })}
      </div>

      {/* HLA + tier */}
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        marginBottom: 8,
        fontFamily: "'Inconsolata', monospace",
        fontSize: 9.5,
      }}>
        <span style={{ color: 'var(--text-meta)' }}>{neo.hla_allele}</span>
        <span style={{ color: tier.color, fontSize: 9 }}>{tier.icon}</span>
      </div>

      {/* Score bar */}
      <div>
        <div style={{
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          marginBottom: 4,
        }}>
          <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 8, color: 'var(--text-muted)' }}>
            composite score
          </span>
          <span style={{ fontFamily: "'Inconsolata', monospace", fontSize: 9, fontWeight: 500, color: scoreColor }}>
            {neo.composite_score.toFixed(3)}
          </span>
        </div>
        <div style={{ height: 3, background: 'var(--border-faint)', borderRadius: 1 }}>
          <div style={{
            height: '100%',
            width: `${neo.composite_score * 100}%`,
            background: scoreColor,
            borderRadius: 1,
          }} />
        </div>
      </div>
    </div>
  );
}

function SmallBadge({ label, color }: { label: string; color: string }) {
  return (
    <span style={{
      fontFamily: "'Inconsolata', monospace",
      fontSize: 8,
      color,
      border: `1px solid ${color}`,
      borderRadius: 2,
      padding: '1px 4px',
      lineHeight: 1,
    }}>
      {label}
    </span>
  );
}

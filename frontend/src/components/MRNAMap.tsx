import type { MRNASegment } from '../api/types';
import { useHover } from '../context/HoverContext';
import { getEpitopeColor } from '../api/helpers';

interface Props {
  segments: MRNASegment[];
  totalLength: number;
}

const STRUCTURAL_COLORS: Record<string, string> = {
  utr5: '#ddd5c8',
  utr3: '#ddd5c8',
  polya: '#c9bfb2',
  linker: 'var(--border-faint)',
};

export default function MRNAMap({ segments, totalLength }: Props) {
  const { hoveredRank, setHoveredRank } = useHover();
  const hasHover = hoveredRank !== null;

  return (
    <div style={{ marginBottom: 20 }}>
      {/* Label */}
      <div style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        letterSpacing: 2.5,
        textTransform: 'uppercase',
        color: 'var(--text-muted)',
        marginBottom: 8,
      }}>
        mRNA Architecture · hover to identify epitope
      </div>

      {/* Track */}
      <div style={{
        display: 'flex',
        height: 28,
        borderRadius: 3,
        overflow: 'hidden',
        border: '1px solid var(--border-light)',
      }}>
        {segments.map((seg, i) => {
          const isEpitope = seg.type === 'epitope';
          const rank = isEpitope && seg.label ? parseInt(seg.label) : null;
          const isHighlighted = rank !== null && hoveredRank === rank;
          const isDimmed = hasHover && !isHighlighted && isEpitope;

          return (
            <div
              key={i}
              style={{
                width: `${seg.pct_of_total}%`,
                background: isEpitope && rank
                  ? getEpitopeColor(rank)
                  : STRUCTURAL_COLORS[seg.type] ?? 'var(--border-faint)',
                opacity: isDimmed ? 0.35 : 1,
                transition: 'opacity 0.15s ease',
                cursor: isEpitope ? 'pointer' : 'default',
                position: 'relative',
              }}
              onMouseEnter={() => rank && setHoveredRank(rank)}
              onMouseLeave={() => rank && setHoveredRank(null)}
              title={isEpitope && seg.label ? `Epitope #${seg.label}` : seg.type}
            />
          );
        })}
      </div>

      {/* Scale line */}
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8,
        color: 'var(--text-faint)',
        marginTop: 4,
      }}>
        <span>1 nt</span>
        <span>{totalLength.toLocaleString()} nt</span>
      </div>

      {/* Legend */}
      <div style={{
        display: 'flex',
        gap: 16,
        marginTop: 8,
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        color: 'var(--text-muted)',
      }}>
        <LegendItem color="#ddd5c8" label="5'/3' UTR" />
        <LegendItem color="#c9bfb2" label="Poly-A" />
        <LegendItem color="var(--border-faint)" label="AAY/GPGPG linkers" />
        <div style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
          <div style={{
            width: 20,
            height: 8,
            borderRadius: 1,
            background: 'linear-gradient(to right, var(--ep-01), var(--ep-10), var(--ep-20))',
          }} />
          Neoantigen epitopes
        </div>
      </div>
    </div>
  );
}

function LegendItem({ color, label }: { color: string; label: string }) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
      <div style={{ width: 8, height: 8, borderRadius: 1, background: color }} />
      {label}
    </div>
  );
}

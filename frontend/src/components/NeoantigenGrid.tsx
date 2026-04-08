import { useState } from 'react';
import type { Neoantigen } from '../api/types';
import NeoantigenCard from './NeoantigenCard';
import NeoantigenModal from './NeoantigenModal';

interface Props {
  neoantigens: Neoantigen[];
}

const INITIAL_SHOW = 9;

export default function NeoantigenGrid({ neoantigens }: Props) {
  const [showAll, setShowAll] = useState(false);
  const [selectedNeo, setSelectedNeo] = useState<Neoantigen | null>(null);
  const visible = showAll ? neoantigens : neoantigens.slice(0, INITIAL_SHOW);
  const hasMore = neoantigens.length > INITIAL_SHOW;

  return (
    <div>
      {/* Legend */}
      <div style={{
        display: 'flex',
        flexWrap: 'wrap',
        gap: 16,
        marginBottom: 8,
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        color: 'var(--text-muted)',
      }}>
        <span><span style={{ color: 'var(--status-complete)' }}>◈</span> AlphaFold</span>
        <span><span style={{ color: '#4a3a1a' }}>◆</span> PANDORA</span>
        <span><span style={{ color: 'var(--text-meta)' }}>◇</span> Position lookup</span>
        <span><span style={{ color: 'var(--status-complete)' }}>fs</span> = frameshift</span>
        <span><span style={{ color: '#1a3a5c' }}>shared</span> = recurrent driver</span>
      </div>

      {/* Peptide key */}
      <div style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        color: 'var(--text-faint)',
        marginBottom: 16,
      }}>
        Peptide key:{' '}
        <span style={{ color: '#8b1a1a', fontWeight: 700, borderBottom: '1.5px solid #8b1a1a' }}>X</span>
        {' '}= mutated residue · {' '}
        <span style={{ color: 'var(--text-secondary)', borderBottom: '1px solid var(--border-light)' }}>X</span>
        {' '}= TCR-facing (P3–P7) · {' '}
        <span style={{ color: 'var(--text-muted)' }}>X</span>
        {' '}= anchor/terminal
      </div>

      {/* Grid */}
      <div style={{
        display: 'grid',
        gridTemplateColumns: 'repeat(3, 1fr)',
        gap: 10,
        alignItems: 'start',
      }}>
        {visible.map(neo => (
          <NeoantigenCard
            key={neo.rank}
            neoantigen={neo}
            onSelect={setSelectedNeo}
          />
        ))}
      </div>

      {/* Show all button */}
      {hasMore && !showAll && (
        <button
          onClick={() => setShowAll(true)}
          style={{
            display: 'block',
            width: '100%',
            marginTop: 12,
            padding: '10px 0',
            fontFamily: "'EB Garamond', serif",
            fontSize: 13,
            color: 'var(--text-meta)',
            background: 'transparent',
            border: '1px solid var(--border-mid)',
            borderRadius: 2,
            cursor: 'pointer',
          }}
        >
          Show all {neoantigens.length} neoantigens ↓
        </button>
      )}

      {/* Modal */}
      {selectedNeo && (
        <NeoantigenModal
          neoantigen={selectedNeo}
          onClose={() => setSelectedNeo(null)}
        />
      )}
    </div>
  );
}

import { useEffect, useRef, useState } from 'react';
// eslint-disable-next-line @typescript-eslint/no-require-imports
import $3Dmol from '3dmol';

interface Props {
  peptide: string;
  hlaAllele: string;
  mutatedPosition: number;
  tcrFacingPositions: number[];
}

// Known pMHC template PDB IDs for common HLA alleles
// These are well-resolved crystal structures we can fetch from RCSB
const HLA_TEMPLATES: Record<string, { pdb: string; chainPeptide: string; chainMHC: string }> = {
  'HLA-A*02:01': { pdb: '1HHI', chainPeptide: 'C', chainMHC: 'A' },
  'HLA-A*03:01': { pdb: '3RL2', chainPeptide: 'C', chainMHC: 'A' },
  'HLA-B*07:02': { pdb: '5EO1', chainPeptide: 'C', chainMHC: 'A' },
  'HLA-B*44:03': { pdb: '1SYS', chainPeptide: 'C', chainMHC: 'A' },
  'HLA-C*05:01': { pdb: '1C1M', chainPeptide: 'P', chainMHC: 'A' },
  'HLA-C*07:02': { pdb: '5VGE', chainPeptide: 'C', chainMHC: 'A' },
};

// Fallback to HLA-A*02:01 (most common, best resolved)
const DEFAULT_TEMPLATE = HLA_TEMPLATES['HLA-A*02:01'];

export default function StructureViewer({ peptide, hlaAllele, mutatedPosition, tcrFacingPositions }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.GLViewer | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const template = HLA_TEMPLATES[hlaAllele] ?? DEFAULT_TEMPLATE;

  useEffect(() => {
    if (!containerRef.current) return;

    setLoading(true);
    setError(null);

    const viewer = $3Dmol.createViewer(containerRef.current, {
      backgroundColor: '#faf7f2',
      antialias: true,
    });
    viewerRef.current = viewer;

    // Fetch PDB from RCSB
    fetch(`https://files.rcsb.org/download/${template.pdb}.pdb`)
      .then(res => {
        if (!res.ok) throw new Error(`PDB fetch failed: ${res.status}`);
        return res.text();
      })
      .then(pdbData => {
        viewer.addModel(pdbData, 'pdb');

        // MHC heavy chain — surface representation, light grey
        viewer.setStyle(
          { chain: template.chainMHC },
          {
            cartoon: { color: '#ddd5c8', opacity: 0.7 },
          }
        );

        // Peptide chain — stick representation
        viewer.setStyle(
          { chain: template.chainPeptide },
          {
            stick: { radius: 0.15, color: '#6b5a47' },
            cartoon: { color: '#9c8b78', opacity: 0.5 },
          }
        );

        // Highlight TCR-facing residues on the peptide in blue
        tcrFacingPositions.forEach(pos => {
          viewer.setStyle(
            { chain: template.chainPeptide, resi: pos + 1 },
            {
              stick: { radius: 0.18, color: '#2e4d8a' },
              cartoon: { color: '#2e4d8a', opacity: 0.6 },
            }
          );
        });

        // Highlight mutated residue in red
        viewer.setStyle(
          { chain: template.chainPeptide, resi: mutatedPosition + 1 },
          {
            stick: { radius: 0.22, color: '#8b1a1a' },
            cartoon: { color: '#8b1a1a' },
          }
        );

        // Add label for mutated residue
        viewer.addLabel(
          peptide[mutatedPosition] + (mutatedPosition + 1),
          {
            position: { x: 0, y: 0, z: 0 },
            fontSize: 11,
            fontColor: '#8b1a1a',
            backgroundColor: 'rgba(255,255,255,0.85)',
            borderColor: '#8b1a1a',
            borderThickness: 1,
            alignment: 'center' as $3Dmol.LabelSpec['alignment'],
          },
          { chain: template.chainPeptide, resi: mutatedPosition + 1 }
        );

        viewer.zoomTo({ chain: template.chainPeptide });
        viewer.zoom(0.7);
        viewer.render();
        setLoading(false);
      })
      .catch(err => {
        setError(err.message);
        setLoading(false);
      });

    return () => {
      viewer.clear();
    };
  }, [peptide, hlaAllele, mutatedPosition, template]);

  return (
    <div style={{
      border: '1px solid var(--border-light)',
      borderRadius: 3,
      overflow: 'hidden',
      marginTop: 10,
    }}>
      {/* Header */}
      <div style={{
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        padding: '6px 10px',
        background: 'var(--bg-sidebar)',
        borderBottom: '1px solid var(--border-light)',
      }}>
        <span style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 10,
          color: 'var(--text-secondary)',
          fontWeight: 500,
        }}>
          pMHC Structure — {hlaAllele}
        </span>
        <span style={{
          fontFamily: "'Inconsolata', monospace",
          fontSize: 8.5,
          color: 'var(--text-muted)',
        }}>
          Template: {template.pdb} (RCSB)
        </span>
      </div>

      {/* Viewer */}
      <div style={{ position: 'relative', height: 300, background: '#faf7f2' }}>
        {loading && (
          <div style={{
            position: 'absolute',
            inset: 0,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            fontFamily: "'Inconsolata', monospace",
            fontSize: 9,
            color: 'var(--text-muted)',
            zIndex: 1,
          }}>
            Loading structure...
          </div>
        )}
        {error && (
          <div style={{
            position: 'absolute',
            inset: 0,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            fontFamily: "'Inconsolata', monospace",
            fontSize: 9,
            color: 'var(--status-error)',
            zIndex: 1,
          }}>
            Could not load structure
          </div>
        )}
        <div
          ref={containerRef}
          style={{ width: '100%', height: '100%' }}
        />
      </div>

      {/* Legend */}
      <div style={{
        display: 'flex',
        gap: 14,
        padding: '6px 10px',
        borderTop: '1px solid var(--border-light)',
        background: 'var(--bg-sidebar)',
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        color: 'var(--text-muted)',
      }}>
        <LegendDot color="#ddd5c8" label="MHC groove" />
        <LegendDot color="#6b5a47" label="Peptide" />
        <LegendDot color="#2e4d8a" label="TCR-facing" />
        <LegendDot color="#8b1a1a" label="Mutated residue" />
      </div>
    </div>
  );
}

function LegendDot({ color, label }: { color: string; label: string }) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
      <div style={{ width: 8, height: 8, borderRadius: '50%', background: color }} />
      {label}
    </div>
  );
}

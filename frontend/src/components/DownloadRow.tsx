import { getDownloadUrl } from '../api/client';

interface Props {
  jobId: string;
}

const DOWNLOAD_FILES = [
  { path: 'module_11/mrna_sequence.fasta', label: 'mrna_sequence.fasta' },
  { path: 'module_11/synthesis_spec.json', label: 'synthesis_spec.json' },
  { path: 'module_09/selected_neoantigens.tsv', label: 'full_report.tsv' },
];

export default function DownloadRow({ jobId }: Props) {
  return (
    <div style={{
      borderTop: '2px solid var(--border-heavy)',
      paddingTop: 20,
      display: 'flex',
      justifyContent: 'space-between',
      alignItems: 'center',
      marginTop: 32,
    }}>
      <span style={{
        fontFamily: "'Spectral', serif",
        fontStyle: 'italic',
        fontSize: 13,
        color: 'var(--text-meta)',
      }}>
        Synthesis-ready outputs — hand to manufacturing facility
      </span>
      <div style={{ display: 'flex', gap: 8 }}>
        {DOWNLOAD_FILES.map(file => (
          <a
            key={file.path}
            href={getDownloadUrl(jobId, file.path)}
            download={file.label}
            style={{
              fontFamily: "'Inconsolata', monospace",
              fontSize: 9.5,
              color: 'var(--text-primary)',
              background: 'transparent',
              border: '1px solid var(--border-mid)',
              padding: '6px 12px',
              borderRadius: 2,
              textDecoration: 'none',
              display: 'inline-flex',
              alignItems: 'center',
              gap: 4,
            }}
          >
            ↓ {file.label}
          </a>
        ))}
      </div>
    </div>
  );
}

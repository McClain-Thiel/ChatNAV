import { useState, useEffect } from 'react';
import { getFilePreview, getDownloadUrl, type FilePreview } from '../api/client';

interface Props {
  jobId: string;
  filePath: string;
  /** Directory prefix for this module's output files */
  outputDir?: string;
}

export default function FileViewer({ jobId, filePath, outputDir }: Props) {
  const [preview, setPreview] = useState<FilePreview | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const fullPath = outputDir ? `${outputDir}/${filePath}` : filePath;

  useEffect(() => {
    setLoading(true);
    setError(null);
    getFilePreview(jobId, fullPath)
      .then(setPreview)
      .catch(err => setError(err.message))
      .finally(() => setLoading(false));
  }, [jobId, fullPath]);

  if (loading) {
    return <div style={styles.container}><span style={styles.meta}>Loading...</span></div>;
  }
  if (error) {
    return <div style={styles.container}><span style={{ ...styles.meta, color: 'var(--status-error)' }}>Could not load file</span></div>;
  }
  if (!preview) return null;

  return (
    <div style={styles.container}>
      {/* File header */}
      <div style={styles.header}>
        <span style={styles.fileName}>{filePath}</span>
        <a
          href={getDownloadUrl(jobId, fullPath)}
          download={filePath}
          style={styles.downloadLink}
        >
          ↓ Download
        </a>
      </div>

      {/* Content */}
      <div style={styles.content}>
        {preview.type === 'table' && <TableView columns={preview.columns!} rows={preview.rows!} />}
        {preview.type === 'json' && <JSONView data={preview.data!} />}
        {preview.type === 'fasta' && <FASTAView sequences={preview.sequences!} />}
        {preview.type === 'text' && <TextView content={preview.content!} />}
      </div>
    </div>
  );
}

function TableView({ columns, rows }: { columns: string[]; rows: string[][] }) {
  const [sortCol, setSortCol] = useState<number | null>(null);
  const [sortAsc, setSortAsc] = useState(true);

  const handleSort = (colIdx: number) => {
    if (sortCol === colIdx) {
      setSortAsc(!sortAsc);
    } else {
      setSortCol(colIdx);
      setSortAsc(true);
    }
  };

  const sortedRows = sortCol !== null
    ? [...rows].sort((a, b) => {
        const va = a[sortCol] ?? '';
        const vb = b[sortCol] ?? '';
        const na = parseFloat(va);
        const nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return sortAsc ? na - nb : nb - na;
        return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      })
    : rows;

  return (
    <div style={{ overflowX: 'auto', maxHeight: 360 }}>
      <table style={styles.table}>
        <thead>
          <tr>
            {columns.map((col, i) => (
              <th
                key={i}
                onClick={() => handleSort(i)}
                style={{
                  ...styles.th,
                  cursor: 'pointer',
                  userSelect: 'none',
                }}
              >
                {col}
                {sortCol === i && (sortAsc ? ' ▲' : ' ▼')}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {sortedRows.map((row, ri) => (
            <tr key={ri} style={{ background: ri % 2 === 0 ? 'transparent' : 'var(--bg-inset)' }}>
              {row.map((cell, ci) => (
                <td key={ci} style={styles.td}>{formatCell(cell)}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function JSONView({ data }: { data: Record<string, unknown> }) {
  // Render as key-value pairs for flat objects, pretty JSON for nested
  const entries = Object.entries(data);
  const isFlat = entries.every(([, v]) => typeof v !== 'object' || v === null);

  if (isFlat) {
    return (
      <div style={{ maxHeight: 360, overflowY: 'auto' }}>
        {entries.map(([key, value]) => (
          <div key={key} style={{
            display: 'flex',
            justifyContent: 'space-between',
            padding: '5px 0',
            borderBottom: '1px solid var(--border-faint)',
          }}>
            <span style={styles.jsonKey}>{key}</span>
            <span style={styles.jsonValue}>{formatValue(value)}</span>
          </div>
        ))}
      </div>
    );
  }

  return (
    <pre style={styles.pre}>
      {JSON.stringify(data, null, 2)}
    </pre>
  );
}

function FASTAView({ sequences }: { sequences: { header: string; sequence: string }[] }) {
  return (
    <div style={{ maxHeight: 360, overflowY: 'auto' }}>
      {sequences.map((seq, i) => (
        <div key={i} style={{ marginBottom: 12 }}>
          {seq.header && (
            <div style={{
              fontFamily: "'Inconsolata', monospace",
              fontSize: 10,
              color: 'var(--text-meta)',
              marginBottom: 4,
            }}>
              &gt;{seq.header}
            </div>
          )}
          <div style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 10.5,
            color: 'var(--text-secondary)',
            lineHeight: 1.8,
            wordBreak: 'break-all',
            letterSpacing: 0.5,
          }}>
            {/* Chunk sequence into lines of 60 characters */}
            {seq.sequence.match(/.{1,60}/g)?.map((chunk, j) => (
              <div key={j}>{chunk}</div>
            ))}
          </div>
        </div>
      ))}
    </div>
  );
}

function TextView({ content }: { content: string }) {
  return <pre style={styles.pre}>{content}</pre>;
}

function formatCell(val: string): string {
  const num = parseFloat(val);
  if (!isNaN(num) && val.includes('.') && num < 100) {
    // Format floating point numbers
    const decimals = val.includes('e') ? 4 : Math.min(val.split('.')[1]?.length ?? 2, 6);
    return num.toFixed(Math.min(decimals, 4));
  }
  // Truncate long strings
  if (val.length > 40) return val.slice(0, 37) + '...';
  return val;
}

function formatValue(val: unknown): string {
  if (val === null || val === undefined) return '—';
  if (typeof val === 'number') {
    if (Number.isInteger(val)) return val.toLocaleString();
    return val.toFixed(4);
  }
  if (typeof val === 'boolean') return val ? 'true' : 'false';
  const s = String(val);
  if (s.length > 60) return s.slice(0, 57) + '...';
  return s;
}

const styles: Record<string, React.CSSProperties> = {
  container: {
    border: '1px solid var(--border-light)',
    borderRadius: 3,
    overflow: 'hidden',
    marginBottom: 12,
  },
  header: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    padding: '6px 10px',
    background: 'var(--bg-sidebar)',
    borderBottom: '1px solid var(--border-light)',
  },
  fileName: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 10,
    color: 'var(--text-secondary)',
    fontWeight: 500,
  },
  downloadLink: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 8.5,
    color: 'var(--text-muted)',
    textDecoration: 'none',
  },
  content: {
    padding: '8px 10px',
    background: 'var(--bg-card)',
  },
  meta: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 9,
    color: 'var(--text-muted)',
    padding: 12,
  },
  table: {
    width: '100%',
    borderCollapse: 'collapse',
    fontFamily: "'Inconsolata', monospace",
    fontSize: 9.5,
  },
  th: {
    textAlign: 'left' as const,
    padding: '4px 8px',
    borderBottom: '1.5px solid var(--border-mid)',
    fontWeight: 500,
    color: 'var(--text-meta)',
    fontSize: 8.5,
    letterSpacing: 0.5,
    textTransform: 'uppercase' as const,
    whiteSpace: 'nowrap' as const,
    position: 'sticky' as const,
    top: 0,
    background: 'var(--bg-card)',
  },
  td: {
    padding: '3px 8px',
    borderBottom: '1px solid var(--border-faint)',
    color: 'var(--text-secondary)',
    whiteSpace: 'nowrap' as const,
  },
  pre: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 10,
    color: 'var(--text-secondary)',
    lineHeight: 1.6,
    maxHeight: 360,
    overflowY: 'auto' as const,
    whiteSpace: 'pre-wrap' as const,
    wordBreak: 'break-all' as const,
    margin: 0,
  },
  jsonKey: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 10,
    color: 'var(--text-meta)',
  },
  jsonValue: {
    fontFamily: "'Inconsolata', monospace",
    fontSize: 10,
    color: 'var(--text-secondary)',
    textAlign: 'right' as const,
    maxWidth: 300,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
  },
};
